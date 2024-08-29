// By: Austin Owens
// Date: 5/23/2024
// Desc: Application to kick off Image FFT Cross Correlation across the system

#include <unistd.h>
#include <fstream>
#include "img_fft_cross_corr.h"
#include <thread>
#include <chrono>

using namespace adf;

int main(int argc, char ** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <xclbin file> <map image file> <iteration>" << std::endl;
        return -1;
    }
    std::cout << "\n##### STARTING APPLICATION #####" << std::endl;

    // User args
    char* xclbin_filename = argv[1];
    char* map_filename = argv[2];
    int iter = std::stoi(argv[3]);

    // Instances of the design instantiated across the PL and AIE
    const int INSTANCES = 4;

    std::cout << "\nLoading xclbin, instantiate PL and graph handles, and malloc memory (HOST)..." << std::endl;
    ImgFFTCrossCorr::startTime();
    ImgFFTCrossCorr ifcc(xclbin_filename, map_filename, iter, INSTANCES);
    ImgFFTCrossCorr::endTime();
    ImgFFTCrossCorr::printTimeDiff("Init completed (HOST)");

    std::cout << "\nReading map image from file (HOST)..." << std::endl;
    ImgFFTCrossCorr::startTime();
    if(!ifcc.openMapFile()) {
        return -1;
    }
    ImgFFTCrossCorr::endTime();
    ImgFFTCrossCorr::printTimeDiff("Read map image completed (HOST)");

    // Run all AIE kernels for specific number of iterations
    std::cout << "\nRun AIE graphs (HOST)... " << std::endl;
    ImgFFTCrossCorr::startTime();
    ifcc.runGraphs();
    ImgFFTCrossCorr::endTime();
    ImgFFTCrossCorr::printTimeDiff("Run AIE graphs completed (HOST)");
    
    // Initialize generic xrt::aie::bo buffer arrays 
    xrt::aie::bo buffers_in[INSTANCES];
    xrt::aie::bo buffers_out[INSTANCES];

    // Initialize matrix multiply xrt::aie::bo buffer arrays 
    xrt::aie::bo buffersA_in[INSTANCES];
    xrt::aie::bo buffersB_in[INSTANCES];

    // Buffer number for telling AIE kernels how many buffers to operating on in parallel
    int buff_num;

    // Structs to hold true and derived peak values
    ImgFFTCrossCorr::PeakValue true_peaks[INSTANCES];
    ImgFFTCrossCorr::PeakValue derived_peaks[INSTANCES];

    std::cout << "\n2D FFT on map image (AIE)..." << std::endl;
    buff_num = 1;
    buffers_in[0] = ifcc.m_map_buffer;
    buffers_out[0] = ifcc.m_map_fft_buffer;
    ImgFFTCrossCorr::startTime();
    ifcc.fft2D(buffers_in, buffers_out, buff_num);
    ImgFFTCrossCorr::endTime();
    ImgFFTCrossCorr::printTimeDiff("2D FFT on map image completed (AIE)");

    for(int n=0; n<iter; n++) {
        std::cout << "\n# ITERATION " << n+1 << ": ";

        std::cout << "\nCreating " << INSTANCES << " template image(s) (HOST)..." << std::endl;
        //ImgFFTCrossCorr::startTime();
        for (int i = 0; i < INSTANCES; i++) {
            true_peaks[i] = ifcc.createTemplateImg(ifcc.m_tmpl_arrays[i]);
        }
        //ImgFFTCrossCorr::endTime();
        //ImgFFTCrossCorr::printTimeDiff("Template image(s) completed (HOST)");

        std::cout << "\n2D FFT on " << INSTANCES << " template image(s) (AIE)..." << std::endl;
        buff_num = INSTANCES;
        for (int i = 0; i < buff_num; i++) {
            buffers_in[i] = ifcc.m_tmpl_buffers[i];
            buffers_out[i] = ifcc.m_tmpl_buffers[i];
        }
        ImgFFTCrossCorr::startTime();
        ifcc.fft2D(buffers_in, buffers_out, buff_num);
        ImgFFTCrossCorr::endTime();
        ImgFFTCrossCorr::printTimeDiff("2D FFT on tmplate image(s) completed (AIE)");

        std::cout << "\nComplex conjugate on " << INSTANCES << " templage image (AIE)..." << std::endl;
        buff_num = INSTANCES;
        for (int i = 0; i < buff_num; i++) {
            buffers_in[i] = ifcc.m_tmpl_buffers[i];
            buffers_out[i] = ifcc.m_tmpl_buffers[i];
        }
        ImgFFTCrossCorr::startTime();
        ifcc.cplxConj(buffers_in, buffers_out, buff_num);
        ImgFFTCrossCorr::endTime();
        ImgFFTCrossCorr::printTimeDiff("Complex conjugate completed (AIE)");

        std::cout << "\nElement-wise matrix multiply between map and template image(s) (AIE)..." << std::endl;
        buff_num = INSTANCES;
        for (int i = 0; i < buff_num; i++) {
            buffersA_in[i] = ifcc.m_map_fft_buffer;
            buffersB_in[i] = ifcc.m_tmpl_buffers[i];
            buffers_out[i] = ifcc.m_tmpl_buffers[i];
        }
        ImgFFTCrossCorr::startTime();
        ifcc.elemMatMult(buffersA_in, buffersB_in, buffers_out, buff_num);
        ImgFFTCrossCorr::endTime();
        ImgFFTCrossCorr::printTimeDiff("Element-wise matrix multiply completed (AIE)");

        std::cout << "\n2D iFFT on " << INSTANCES << " image(s) (AIE)..." << std::endl;
        buff_num = INSTANCES;
        for (int i = 0; i < buff_num; i++) {
            buffers_in[i] = ifcc.m_tmpl_buffers[i];
            buffers_out[i] = ifcc.m_tmpl_buffers[i];
        }
        ImgFFTCrossCorr::startTime();
        ifcc.ifft2D(buffers_in, buffers_out, buff_num);
        ImgFFTCrossCorr::endTime();
        ImgFFTCrossCorr::printTimeDiff("2D iFFT on image(s) completed (AIE)");

        std::cout << "\nPeak search on " << INSTANCES << " image(s) (AIE)..." << std::endl;
        buff_num = INSTANCES;
        TT_DATA tmp_prev_peak_vals[buff_num];
        for (int i = 0; i < buff_num; i++) {
            tmp_prev_peak_vals[i] = ifcc.m_peak_arrays[i][0];
            buffers_in[i] = ifcc.m_tmpl_buffers[i];
            buffers_out[i] = ifcc.m_peak_buffers[i];
        }
        ImgFFTCrossCorr::startTime();
        ifcc.peakSearch(buffers_in, buffers_out, buff_num);
        // TODO: This is a temporary fix. I should find out why wait() blocking operation
        // inside ifcc.peakSearch is not working
        
        for (int i = 0; i < buff_num; i++) {
            while (tmp_prev_peak_vals[i].imag == ifcc.m_peak_arrays[i][0].imag) {
                std::this_thread::sleep_for(std::chrono::milliseconds(1));
                //std::cout << "prev imag: " << ifcc.m_peak_arrays[i][0].imag << "| imag: " << tmp_prev_peak_vals[i].imag <<std::endl;
            }
        }
        ImgFFTCrossCorr::endTime();
        ImgFFTCrossCorr::printTimeDiff("Peak search on image(s) completed (AIE)");

        std::cout << "\nPost peak search offset on " << INSTANCES << " images(s) (HOST)..." << std::endl;
        ImgFFTCrossCorr::startTime();
        for (int i = 0; i < INSTANCES; i++) {
            derived_peaks[i] = ifcc.postPeakSearchOffset(ifcc.m_peak_arrays[i][0]);
        }
        ImgFFTCrossCorr::endTime();
        ImgFFTCrossCorr::printTimeDiff("Post peak search offset on images(s) completed (HOST)");
        
        for (int i = 0; i < INSTANCES; i++) {
            ifcc.printResults(true_peaks[i], derived_peaks[i]);
        }
        ImgFFTCrossCorr::printTotalTime(n);
    }

    ImgFFTCrossCorr::printAvgTime(iter);
    
    ifcc.endGraphs();

    return 0;
};
