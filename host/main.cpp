// By: Austin Owens
// Date: 5/23/2024
// Desc: Application to kick off Image FFT Cross Correlation across the system

#include <unistd.h>
#include <fstream>
#include "sar_backproject.h"
#include <thread>
#include <chrono>

using namespace adf;

//void print_arr(TT_DATA* arr, int num_rows, int num_cols, int divisor=1) {
//    for(int r=0; r<num_rows; r++) {
//        for(int c=0; c<num_cols; c++) {
//            int index = (r*MAT_COLS) + c;
//            printf("array[%d][%d] = {%f, %f}\n", r, c, arr[index].real/divisor, arr[index].imag/divisor);
//        }
//    }
//}
//
//int ifftErrorCheck(TT_DATA* data_array) {
//    int error_cnt = 0;
//    float real, imag;
//
//    for(int r = 0; r < MAT_ROWS; r++) {
//        for(int c = 0; c < MAT_COLS; c++) {
//            int index = (r*MAT_COLS)+c;
//            real = data_array[index].real/TP_POINT_SIZE;
//            imag = data_array[index].imag/TP_POINT_SIZE;
//            //printf("IFFT[%d, %d] = {%f, %f}\n", r, c, real, imag);
//
//            if (c == 0) {
//                if(std::abs(real-1.5) > 0.00001 || std::abs(imag+1.5) > 0.00001) {
//                    printf("ERROR IFFT[%d, %d] = {%lf, %lf} | EXPECTED {%f, %f}\n", r, c, real, imag, 1.5, 1.5);
//                    error_cnt++;
//                }
//            }
//            else {
//                if(std::abs(real) > 0.00001 || std::abs(imag) > 0.00001) {
//                    printf("ERROR IFFT[%d, %d] = {%lf, %lf} | EXPECTED {%d, %d}\n", r, c, real, imag, 0, 0);
//                    error_cnt++;
//                }
//
//            }
//        }
//    }
//
//    printf("IFFT ERRORS: %d\n\n", error_cnt);
//
//    return error_cnt;
//}

int main(int argc, char ** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <xclbin file> <slowtime dataset csv file> <range compressed dataset csv file> <img out csv file> <iteration>" << std::endl;
        return -1;
    }
    std::cout << "\n##### STARTING APPLICATION #####" << std::endl;

    // User args
    char* xclbin_filename = argv[1];
    char* st_dataset_filename = argv[2];
    char* rc_dataset_filename = argv[3];
    char* img_out_filename = argv[4];
    int iter = std::stoi(argv[5]);

    // Instances of the design instantiated across the PL and AIE
    const int INSTANCES = 1;

    std::cout << "\nLoading xclbin, instantiate AIE graph handles, and malloc memory (HOST)..." << std::endl;
    SARBackproject::startTime();
    SARBackproject ifcc(xclbin_filename, st_dataset_filename, rc_dataset_filename, img_out_filename, iter, INSTANCES);
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Init completed (HOST)");

    std::cout << "\nPopulating data buffers (HOST)..." << std::endl;
    SARBackproject::startTime();
    if(ifcc.fetchRadarData() != 0) {
        std::cout << "\nPopulating data buffers failed (HOST)" << std::endl;
        return -1;
    }
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Populating data buffers completed (HOST)");

    std::cout << "\nGenerating target pixels (HOST)..." << std::endl;
    SARBackproject::startTime();
    ifcc.genTargetPixels();
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Generating target pixels  completed (HOST)");

    // Start all AIE kernels
    std::cout << "\nRun AIE graphs (HOST)... " << std::endl;
    SARBackproject::startTime();
    ifcc.runGraphs();
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Run AIE graphs completed (HOST)");

    // Initialize bp xrt::aie::bo buffer arrays 
    xrt::aie::bo buffers_broadcast_data_in[INSTANCES];
    //xrt::aie::bo buffers_x_ant_pos_in[INSTANCES];
    //xrt::aie::bo buffers_y_ant_pos_in[INSTANCES];
    //xrt::aie::bo buffers_z_ant_pos_in[INSTANCES];
    //xrt::aie::bo buffers_ref_range_in[INSTANCES];
    xrt::aie::bo buffers_rc_in[INSTANCES];
    xrt::aie::bo buffers_xyz_px_in[INSTANCES];
    //xrt::aie::bo buffers_img_out[INSTANCES];
    xrt::bo buffers_img_out[INSTANCES];

    // Initialize generic xrt::aie::bo buffer arrays 
    //xrt::aie::bo buffers_in[INSTANCES];
    //xrt::aie::bo buffers_out[INSTANCES];

    // Initialize matrix multiply xrt::aie::bo buffer arrays 
    //xrt::aie::bo buffersA_in[INSTANCES];
    //xrt::aie::bo buffersB_in[INSTANCES];

    // Buffer number for telling AIE kernels how many buffers to operate on in parallel
    int buff_num = 1;

    std::cout << "\nPerform Backprojection (AIE)..." << std::endl;
    buffers_broadcast_data_in[0] = ifcc.m_broadcast_data_buffer;
    //buffers_x_ant_pos_in[0] = ifcc.m_x_ant_pos_buffer;
    //buffers_y_ant_pos_in[0] = ifcc.m_y_ant_pos_buffer;
    //buffers_z_ant_pos_in[0] = ifcc.m_z_ant_pos_buffer;
    //buffers_ref_range_in[0] = ifcc.m_ref_range_buffer;
    buffers_rc_in[0] = ifcc.m_rc_buffer;
    buffers_xyz_px_in[0] = ifcc.m_xyz_px_buffer;
    //buffers_img_out[0] = ifcc.m_img_buffer;
    buffers_img_out[0] = ifcc.m_dma_pkt_router_buffers[0];
    SARBackproject::startTime();
    ifcc.bp(buffers_broadcast_data_in, buffers_rc_in, 
            buffers_xyz_px_in, buffers_img_out, buff_num);
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Backprojection completed (AIE)");
    
    if(ifcc.writeImg() != 0) {
        std::cout << "\nWriting image failed!" << std::endl;
        return -1;
    }

    //for(int i=0; i<64; i++) {
    //    printf("%d: %f\t%f\n", i, ifcc.m_img_array[i].real, ifcc.m_img_array[i].imag);
    //}
    //for(int i=0; i<4; i++) {
    //    for(int j=i*2048; j<(i*2048)+8; j++) {
    //        printf("%d: %f\t%f\n", j, ifcc.m_img_array[j].real, ifcc.m_img_array[j].imag);
    //    }
    //}

    //std::cout << "\nReshape reference function for FFT (HOST - OpenMP)..." << std::endl;
    //SARBackproject::startTime();
    //ifcc.reshapeMatrix(ifcc.m_ref_func_array, MAT_ROWS, MAT_COLS, FFT_NPORTS);
    //SARBackproject::endTime();
    //SARBackproject::printTimeDiff("Reshape completed (HOST - OpenMP)");

    //std::cout << "\nFFT on reference function (AIE)..." << std::endl;
    //buffers_in[0] = ifcc.m_ref_func_buffer;
    //buffers_out[0] = ifcc.m_ref_func_buffer;
    //SARBackproject::startTime();
    //ifcc.fft(buffers_in, buffers_out, buff_num);
    //SARBackproject::endTime();
    //SARBackproject::printTimeDiff("FFT completed (AIE)");
    //
    ////TODO: DEBUG
    ////print_arr(ifcc.m_ref_func_array, 5, 5);

    //std::cout << "\nComplex conjugate on reference function (AIE)..." << std::endl;
    //buffers_in[0] = ifcc.m_ref_func_buffer;
    //buffers_out[0] = ifcc.m_ref_func_buffer;
    //SARBackproject::startTime();
    //ifcc.cplxConj(buffers_in, buffers_out, buff_num);
    //SARBackproject::endTime();
    //SARBackproject::printTimeDiff("Complex conjugate completed (AIE)");

    ////TODO: DEBUG
    ////print_arr(ifcc.m_ref_func_array, 5, 5);

    //std::cout << "\nReshape range data for FFT (HOST - OpenMP)..." << std::endl;
    //SARBackproject::startTime();
    //ifcc.reshapeMatrix(ifcc.m_range_data_array, MAT_ROWS, MAT_COLS, FFT_NPORTS);
    //SARBackproject::endTime();
    //SARBackproject::printTimeDiff("Reshape completed (HOST - OpenMP)");

    //std::cout << "\nFFT on range data (AIE)..." << std::endl;
    //buffers_in[0] = ifcc.m_range_data_buffer;
    //buffers_out[0] = ifcc.m_range_data_buffer;
    //SARBackproject::startTime();
    //ifcc.fft(buffers_in, buffers_out, buff_num);
    //SARBackproject::endTime();
    //SARBackproject::printTimeDiff("FFT completed (AIE)");

    ////TODO: DEBUG
    ////print_arr(ifcc.m_range_data_array, 5, 5);

    //std::cout << "\nElement-wise matrix multiply between ref function and range data (AIE)..." << std::endl;
    //buffersA_in[0] = ifcc.m_range_data_buffer;
    //buffersB_in[0] = ifcc.m_ref_func_buffer;
    //buffers_out[0] = ifcc.m_range_data_buffer;
    //SARBackproject::startTime();
    //ifcc.elemMatMult(buffersA_in, buffersB_in, buffers_out, buff_num);
    //SARBackproject::endTime();
    //SARBackproject::printTimeDiff("Element-wise matrix multiply completed (AIE)");

    ////TODO: DEBUG
    ////print_arr(ifcc.m_range_data_array, 5, 5);

    //std::cout << "\niFFT on range data (AIE)..." << std::endl;
    //buffers_in[0] = ifcc.m_range_data_buffer;
    //buffers_out[0] = ifcc.m_range_data_buffer;
    //SARBackproject::startTime();
    //ifcc.ifft(buffers_in, buffers_out, buff_num);
    //SARBackproject::endTime();
    //SARBackproject::printTimeDiff("iFFT completed (AIE)");

    ////TODO: DEBUG
    ////print_arr(ifcc.m_range_data_array, 5, 5, TP_POINT_SIZE);
    ////ifftErrorCheck(ifcc.m_range_data_array);

    //std::cout << "\nReshape range data after IFFT (HOST - OpenMP)..." << std::endl;
    //SARBackproject::startTime();
    //ifcc.reshapeMatrix(ifcc.m_range_data_array, MAT_ROWS, MAT_COLS, FFT_NPORTS, true);
    //SARBackproject::endTime();
    //SARBackproject::printTimeDiff("Reshape completed (HOST - OpenMP)");

    //TT_DATA bp_ref;
    //for(int az=0; az<MAT_ROWS; az++) {
    //    bp_ref = exp((j*4*pi*az)/lambda);

    //    std::cout << "\nFFT on reference function (AIE)..." << std::endl;
    //    buffers_in[0] = ifcc.m_bp_ref_func_buffer;
    //    buffers_out[0] = ifcc.m_bp_ref_func_buffer;
    //    SARBackproject::startTime();
    //    ifcc.fft(buffers_in, buffers_out, buff_num);
    //    SARBackproject::endTime();
    //    SARBackproject::printTimeDiff("FFT completed (AIE)");

    //}

    

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

    return 0;
};
