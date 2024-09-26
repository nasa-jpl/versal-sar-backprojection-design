// By: Austin Owens
// Date: 5/23/2024
// Desc: Application to kick off Image FFT Cross Correlation across the system

#include <unistd.h>
#include <fstream>
#include "sar_backproject.h"
#include <thread>
#include <chrono>

using namespace adf;

int main(int argc, char ** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <xclbin file> <hdf5 file> <iteration>" << std::endl;
        return -1;
    }
    std::cout << "\n##### STARTING APPLICATION #####" << std::endl;

    // User args
    char* xclbin_filename = argv[1];
    char* hdf5_filename = argv[2];
    int iter = std::stoi(argv[3]);

    // Instances of the design instantiated across the PL and AIE
    const int INSTANCES = 1;

    std::cout << "\nLoading xclbin, instantiate AIE graph handles, and malloc memory (HOST)..." << std::endl;
    SARBackproject::startTime();
    SARBackproject ifcc(xclbin_filename, hdf5_filename, iter, INSTANCES);
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Init completed (HOST)");

    std::cout << "\nReading range data from file (HOST)..." << std::endl;
    SARBackproject::startTime();
    if(!ifcc.openHDF5File()) {
        return -1;
    }
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Read range data completed (HOST)");

    // Run all AIE kernels for specific number of iterations
    std::cout << "\nRun AIE graphs (HOST)... " << std::endl;
    SARBackproject::startTime();
    ifcc.runGraphs();
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Run AIE graphs completed (HOST)");
    
    // Initialize generic xrt::aie::bo buffer arrays 
    xrt::aie::bo buffers_in[INSTANCES];
    xrt::aie::bo buffers_out[INSTANCES];

    // Initialize matrix multiply xrt::aie::bo buffer arrays 
    xrt::aie::bo buffersA_in[INSTANCES];
    xrt::aie::bo buffersB_in[INSTANCES];

    // Buffer number for telling AIE kernels how many buffers to operate on in parallel
    int buff_num;

    std::cout << "\nFFT on range data (AIE)..." << std::endl;
    buff_num = 1;
    buffers_in[0] = ifcc.m_range_data_buffer;
    buffers_out[0] = ifcc.m_range_data_buffer;
    SARBackproject::startTime();
    ifcc.fft(buffers_in, buffers_out, buff_num);
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("FFT on range data completed (AIE)");

    //TODO: TEMPORARY
    for(int r=0; r<5; r++) {
        for(int c=0; c<5; c++) {
            int index = (r*MAT_COLS) + c;
            printf("ifcc.m_range_data_array[%d][%d] = {%f, %f}\n", r, c, ifcc.m_range_data_array[index].real, ifcc.m_range_data_array[index].imag);
        }
    }

    //std::cout << "\nComplex conjugate on range data (AIE)..." << std::endl;
    //buff_num = INSTANCES;
    //for (int i = 0; i < buff_num; i++) {
    //    buffers_in[i] = ifcc.m_tmpl_buffers[i];
    //    buffers_out[i] = ifcc.m_tmpl_buffers[i];
    //}
    //SARBackproject::startTime();
    //ifcc.cplxConj(buffers_in, buffers_out, buff_num);
    //SARBackproject::endTime();
    //SARBackproject::printTimeDiff("Complex conjugate completed (AIE)");

    std::cout << "\nElement-wise matrix multiply (AIE)..." << std::endl;
    buffersA_in[0] = ifcc.m_ref_func_range_comp_buffer;
    buffersB_in[0] = ifcc.m_range_data_buffer;
    buffers_out[0] = ifcc.m_range_data_buffer;
    SARBackproject::startTime();
    ifcc.elemMatMult(buffersA_in, buffersB_in, buffers_out, buff_num);
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Element-wise matrix multiply completed (AIE)");

    //TODO: TEMPORARY
    for(int r=0; r<5; r++) {
        for(int c=0; c<5; c++) {
            int index = (r*MAT_COLS) + c;
            printf("ifcc.m_range_data_array[%d][%d] = {%f, %f}\n", r, c, ifcc.m_range_data_array[index].real, ifcc.m_range_data_array[index].imag);
        }
    }

    std::cout << "\niFFT on range data (AIE)..." << std::endl;
    buffers_in[0] = ifcc.m_range_data_buffer;
    buffers_out[0] = ifcc.m_range_data_buffer;
    SARBackproject::startTime();
    ifcc.ifft(buffers_in, buffers_out, buff_num);
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("iFFT on range data completed (AIE)");

    //TODO: TEMPORARY
    for(int r=0; r<5; r++) {
        for(int c=0; c<5; c++) {
            int index = (r*MAT_COLS) + c;
            printf("ifcc.m_range_data_array[%d][%d] = {%f, %f}\n", r, c, ifcc.m_range_data_array[index].real/TP_POINT_SIZE, ifcc.m_range_data_array[index].imag/TP_POINT_SIZE);
        }
    }

    //for(int n=0; n<iter; n++) {
    //    std::cout << "\n# ITERATION " << n+1 << ": ";

    //    std::cout << "\nCreating " << INSTANCES << " template image(s) (HOST)..." << std::endl;
    //    //SARBackproject::startTime();
    //    for (int i = 0; i < INSTANCES; i++) {
    //        true_peaks[i] = ifcc.createTemplateImg(ifcc.m_tmpl_arrays[i]);
    //    }
    //    //SARBackproject::endTime();
    //    //SARBackproject::printTimeDiff("Template image(s) completed (HOST)");

    //    std::cout << "\n2D FFT on " << INSTANCES << " template image(s) (AIE)..." << std::endl;
    //    buff_num = INSTANCES;
    //    for (int i = 0; i < buff_num; i++) {
    //        buffers_in[i] = ifcc.m_tmpl_buffers[i];
    //        buffers_out[i] = ifcc.m_tmpl_buffers[i];
    //    }
    //    SARBackproject::startTime();
    //    ifcc.fft2D(buffers_in, buffers_out, buff_num);
    //    SARBackproject::endTime();
    //    SARBackproject::printTimeDiff("2D FFT on tmplate image(s) completed (AIE)");

    //    std::cout << "\nComplex conjugate on " << INSTANCES << " templage image (AIE)..." << std::endl;
    //    buff_num = INSTANCES;
    //    for (int i = 0; i < buff_num; i++) {
    //        buffers_in[i] = ifcc.m_tmpl_buffers[i];
    //        buffers_out[i] = ifcc.m_tmpl_buffers[i];
    //    }
    //    SARBackproject::startTime();
    //    ifcc.cplxConj(buffers_in, buffers_out, buff_num);
    //    SARBackproject::endTime();
    //    SARBackproject::printTimeDiff("Complex conjugate completed (AIE)");

    //    std::cout << "\nElement-wise matrix multiply between map and template image(s) (AIE)..." << std::endl;
    //    buff_num = INSTANCES;
    //    for (int i = 0; i < buff_num; i++) {
    //        buffersA_in[i] = ifcc.m_map_fft_buffer;
    //        buffersB_in[i] = ifcc.m_tmpl_buffers[i];
    //        buffers_out[i] = ifcc.m_tmpl_buffers[i];
    //    }
    //    SARBackproject::startTime();
    //    ifcc.elemMatMult(buffersA_in, buffersB_in, buffers_out, buff_num);
    //    SARBackproject::endTime();
    //    SARBackproject::printTimeDiff("Element-wise matrix multiply completed (AIE)");

    //    std::cout << "\n2D iFFT on " << INSTANCES << " image(s) (AIE)..." << std::endl;
    //    buff_num = INSTANCES;
    //    for (int i = 0; i < buff_num; i++) {
    //        buffers_in[i] = ifcc.m_tmpl_buffers[i];
    //        buffers_out[i] = ifcc.m_tmpl_buffers[i];
    //    }
    //    SARBackproject::startTime();
    //    ifcc.ifft2D(buffers_in, buffers_out, buff_num);
    //    SARBackproject::endTime();
    //    SARBackproject::printTimeDiff("2D iFFT on image(s) completed (AIE)");

    //    std::cout << "\nPeak search on " << INSTANCES << " image(s) (AIE)..." << std::endl;
    //    buff_num = INSTANCES;
    //    TT_DATA tmp_prev_peak_vals[buff_num];
    //    for (int i = 0; i < buff_num; i++) {
    //        tmp_prev_peak_vals[i] = ifcc.m_peak_arrays[i][0];
    //        buffers_in[i] = ifcc.m_tmpl_buffers[i];
    //        buffers_out[i] = ifcc.m_peak_buffers[i];
    //    }
    //    SARBackproject::startTime();
    //    ifcc.peakSearch(buffers_in, buffers_out, buff_num);
    //    // TODO: This is a temporary fix. I should find out why wait() blocking operation
    //    // inside ifcc.peakSearch is not working
    //    
    //    for (int i = 0; i < buff_num; i++) {
    //        while (tmp_prev_peak_vals[i].imag == ifcc.m_peak_arrays[i][0].imag) {
    //            std::this_thread::sleep_for(std::chrono::milliseconds(1));
    //            //std::cout << "prev imag: " << ifcc.m_peak_arrays[i][0].imag << "| imag: " << tmp_prev_peak_vals[i].imag <<std::endl;
    //        }
    //    }
    //    SARBackproject::endTime();
    //    SARBackproject::printTimeDiff("Peak search on image(s) completed (AIE)");

    //    std::cout << "\nPost peak search offset on " << INSTANCES << " images(s) (HOST)..." << std::endl;
    //    SARBackproject::startTime();
    //    for (int i = 0; i < INSTANCES; i++) {
    //        derived_peaks[i] = ifcc.postPeakSearchOffset(ifcc.m_peak_arrays[i][0]);
    //    }
    //    SARBackproject::endTime();
    //    SARBackproject::printTimeDiff("Post peak search offset on images(s) completed (HOST)");
    //    
    //    for (int i = 0; i < INSTANCES; i++) {
    //        ifcc.printResults(true_peaks[i], derived_peaks[i]);
    //    }
    //    SARBackproject::printTotalTime(n);
    //}

    //SARBackproject::printAvgTime(iter);
    
    ifcc.endGraphs();

    return 0;
};
