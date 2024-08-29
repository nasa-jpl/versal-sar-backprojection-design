// By: Austin Owens
// Date: 6/27/2024
// Desc: Class to handle image FFT cross correlation 

#include "img_fft_cross_corr.h"
#include <omp.h>
#include <fstream>
#include "../common.h"
#include <thread>
#include <chrono>

double ImgFFTCrossCorr::total_time = 0;
double ImgFFTCrossCorr::total_avg_time = 0;
struct timespec ImgFFTCrossCorr::time_start = {0, 0};
struct timespec ImgFFTCrossCorr::time_end = {0, 0};

ImgFFTCrossCorr::ImgFFTCrossCorr(const char* xclbin_filename, const char* map_filename, int iter, int instances)
: m_xclbin_filename(xclbin_filename)
, m_map_filename(map_filename)
, m_iter(iter)
, m_instances(instances)
, m_device(0)
, m_uuid(m_device.load_xclbin(this->m_xclbin_filename))
, m_map_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_map_array(m_map_buffer.map<TT_DATA*>())
, m_map_fft_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_map_fft_array(m_map_fft_buffer.map<TT_DATA*>())
{
    // Instantiate kernels, handlers, and buffers for multiple instances
    for(int i = 0; i < this->m_instances; i++) {
        std::string dma_hls_fft_kernel_str = "dma_hls:{dma_hls_" + std::to_string(i * 2) + "}";
        std::string dma_hls_ifft_kernel_str = "dma_hls:{dma_hls_" + std::to_string(i * 2 + 1) + "}";
        std::string fftRowsGraph_str = "fftRowsGraph[" + std::to_string(i) + "]";
        std::string fftColsGraph_str = "fftColsGraph[" + std::to_string(i) + "]";
        std::string ifftColsGraph_str = "ifftColsGraph[" + std::to_string(i) + "]";
        std::string ifftRowsGraph_str = "ifftRowsGraph[" + std::to_string(i) + "]";
        std::string cplxConjGraph_str = "cplxConjGraph[" + std::to_string(i) + "]";
        std::string hpGraph_str = "hpGraph[" + std::to_string(i) + "]";
        std::string peakGraph_str = "peakGraph[" + std::to_string(i) + "]";

        m_dma_hls_fft_kernels.push_back(xrt::kernel(m_device, m_uuid, dma_hls_fft_kernel_str));
        m_dma_hls_fft_run_hdls.push_back(xrt::run(m_dma_hls_fft_kernels[i]));
        m_dma_hls_fft_buffers.push_back(xrt::bo(m_device, BLOCK_SIZE_BYTES, m_dma_hls_fft_kernels[i].group_id(0)));
        m_dma_hls_ifft_kernels.push_back(xrt::kernel(m_device, m_uuid, dma_hls_ifft_kernel_str));
        m_dma_hls_ifft_run_hdls.push_back(xrt::run(m_dma_hls_ifft_kernels[i]));
        m_dma_hls_ifft_buffers.push_back(xrt::bo(m_device, BLOCK_SIZE_BYTES, m_dma_hls_ifft_kernels[i].group_id(0)));
        m_fft_rows_graph_hdls.push_back(xrt::graph(m_device, m_uuid, fftRowsGraph_str));
        m_fft_cols_graph_hdls.push_back(xrt::graph(m_device, m_uuid, fftColsGraph_str));
        m_ifft_cols_graph_hdls.push_back(xrt::graph(m_device, m_uuid, ifftColsGraph_str));
        m_ifft_rows_graph_hdls.push_back(xrt::graph(m_device, m_uuid, ifftRowsGraph_str));
        m_cplx_conj_graph_hdls.push_back(xrt::graph(m_device, m_uuid, cplxConjGraph_str));
        m_hp_graph_hdls.push_back(xrt::graph(m_device, m_uuid, hpGraph_str));
        m_peak_graph_hdls.push_back(xrt::graph(m_device, m_uuid, peakGraph_str));
        m_aie_to_pl_buffers.push_back(xrt::aie::bo(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0));
        m_tmpl_buffers.push_back(xrt::aie::bo(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0));
        m_tmpl_arrays.push_back(m_tmpl_buffers[i].map<TT_DATA*>());
        m_peak_buffers.push_back(xrt::aie::bo(m_device, sizeof(TT_DATA), xrt::bo::flags::normal, 0));
        m_peak_arrays.push_back(m_peak_buffers[i].map<TT_DATA*>());
    }

    // Seed random number generator used for selecting template images in map image
    std::time_t seed = std::time(0);
    //std::time_t seed = 1719629552;
    std::srand(seed);
    std::cout << "Rand Seed: " << seed << std::endl;
}

void ImgFFTCrossCorr::startTime() {
    clock_gettime(CLOCK_MONOTONIC, &time_start);
}

void ImgFFTCrossCorr::endTime() {
    clock_gettime(CLOCK_MONOTONIC, &time_end);
}

void ImgFFTCrossCorr::printTimeDiff(const char *msg) {
    double elapsed_time_millis;

    // Calculating the elapsed time in milliseconds

    // Seconds to milliseconds
    elapsed_time_millis = (time_end.tv_sec - time_start.tv_sec) * 1000.0;

    // Nanoseconds to milliseconds
    elapsed_time_millis += (time_end.tv_nsec - time_start.tv_nsec) / 1000000.0;

    // Save accumulated time
    total_time += elapsed_time_millis;

    printf("Elapsed time (%s): %.2f milliseconds\n", msg, elapsed_time_millis);
}

void ImgFFTCrossCorr::printTotalTime(int curr_iter) {
    printf("Elapsed total time: %.2f milliseconds\n", total_time);

    // Don't include the first iteration since it includes calculating the map image
    if(curr_iter > 0) {
        total_avg_time += total_time;
    }
    total_time = 0;
}

void ImgFFTCrossCorr::printAvgTime(int iterations) {
    if (iterations > 1) {
        total_avg_time = total_avg_time/(iterations-1);
        printf("Elapsed total average time over this last %d iteration(s): %.2f milliseconds\n", iterations-1, total_avg_time);
    }
}

void ImgFFTCrossCorr::runGraphs() {
    this->m_fft_rows_graph_hdls[0].run(MAT_ROWS + (MAT_ROWS*this->m_iter));
    this->m_fft_cols_graph_hdls[0].run(MAT_COLS + (MAT_COLS*this->m_iter));
    this->m_ifft_cols_graph_hdls[0].run(MAT_COLS*this->m_iter);
    this->m_ifft_rows_graph_hdls[0].run(MAT_ROWS*this->m_iter);
    this->m_cplx_conj_graph_hdls[0].run(MAT_ROWS*this->m_iter);
    this->m_hp_graph_hdls[0].run((TP_DIM/TP_NUM_FRAMES)*this->m_iter);
    this->m_peak_graph_hdls[0].run(this->m_iter);

    for(int i = 1; i < this->m_instances; i++) {
        this->m_fft_rows_graph_hdls[i].run(MAT_ROWS*this->m_iter);
        this->m_fft_cols_graph_hdls[i].run(MAT_COLS*this->m_iter);
        this->m_ifft_cols_graph_hdls[i].run(MAT_COLS*this->m_iter);
        this->m_ifft_rows_graph_hdls[i].run(MAT_ROWS*this->m_iter);
        this->m_cplx_conj_graph_hdls[i].run(MAT_ROWS*this->m_iter);
        this->m_hp_graph_hdls[i].run((TP_DIM/TP_NUM_FRAMES)*this->m_iter);
        this->m_peak_graph_hdls[i].run(this->m_iter);
    }
}

void ImgFFTCrossCorr::endGraphs() {
    for(int i = 0; i < this->m_instances; i++) {
        this->m_fft_rows_graph_hdls[i].end();
        this->m_fft_cols_graph_hdls[i].end();
        this->m_ifft_cols_graph_hdls[i].end();
        this->m_ifft_rows_graph_hdls[i].end();
        this->m_cplx_conj_graph_hdls[i].end();
        this->m_hp_graph_hdls[i].end();
        this->m_peak_graph_hdls[i].end();
    }
}

bool ImgFFTCrossCorr::openMapFile() {
    std::ifstream map_file(this->m_map_filename);
    if (!map_file.is_open()) {
        std::cerr << "Error opening map file: " << this->m_map_filename << std::endl;
        return false;
    }

    for (int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
        std::string hex_str;
        if (!(map_file >> hex_str)) {
            std::cerr << "Error reading map image file: insufficient data" << std::endl;
            map_file.close();
            return false;
        }
        int val = std::stoi(hex_str, nullptr, 16);
        if (val >= 0x8000) {
            val |= 0xFFFF0000;
        }

        float real = static_cast<float>(val + 0x80);
        this->m_map_array[i] = (TT_DATA) {real, 0};
    }

    map_file.close();

    return true;
}

ImgFFTCrossCorr::PeakValue ImgFFTCrossCorr::createTemplateImg(TT_DATA* array_in) {
    PeakValue true_peak;

    for (int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
        array_in[i] = (TT_DATA) {0, 0};
    }

    // Calculate the top-left corner location for the template region
    int map_start_row = std::rand() % (MAT_ROWS - TMPL_MAT_ROWS + 1);
    int map_start_col = std::rand() % (MAT_COLS - TMPL_MAT_COLS + 1);

    true_peak.x = map_start_row + (TMPL_MAT_ROWS/2);
    true_peak.y = map_start_col + (TMPL_MAT_COLS/2);
    printf("True Peak Info: peak.x = %d | peak.y = %d\n", true_peak.x, true_peak.y);

    // Calculate the top-left corner of what will be the template image 
    // (which will always be in the center)
    int tmpl_start_row = (MAT_ROWS - TMPL_MAT_ROWS) / 2;
    int tmpl_start_col = (MAT_COLS - TMPL_MAT_COLS) / 2;

    // Copy smaller region from map image into center of template image
    int map_index, tmpl_index;
    for (int i = 0; i < TMPL_MAT_ROWS; ++i) {
        for (int j = 0; j < TMPL_MAT_COLS; ++j) {
            map_index = (map_start_row + i) * MAT_COLS + (map_start_col + j);
            tmpl_index = (tmpl_start_row + i) * MAT_COLS + (tmpl_start_col + j);
            array_in[tmpl_index] = this->m_map_array[map_index];
        }
    }

    return true_peak;
}

void ImgFFTCrossCorr::fft2D(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
    int offset = 0;
    std::vector<xrt::bo::async_handle> rows_buff_async_hdls;
    std::vector<xrt::bo::async_handle> cols_buff_async_hdls;

    for(int i = 0; i < num_of_buffers; i++) {
        buffers_in[i].async("fftRowsGraph[" + std::to_string(i) + "].gmio_in", 
                            XCL_BO_SYNC_BO_GMIO_TO_AIE, BLOCK_SIZE_BYTES, offset);

        rows_buff_async_hdls.push_back(m_aie_to_pl_buffers[i].async("fftRowsGraph[" + std::to_string(i) + "].gmio_out", 
                                                            XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                            BLOCK_SIZE_BYTES, 
                                                            offset));
    }


    // Setup and run PL kernels for data transfer to AIE
    for(int i = 0; i < num_of_buffers; i++) {
        rows_buff_async_hdls[i].wait();
        this->m_dma_hls_fft_buffers[i] = xrt::bo(m_aie_to_pl_buffers[i]);
        this->m_dma_hls_fft_buffers[i].sync(XCL_BO_SYNC_BO_TO_DEVICE, BLOCK_SIZE_BYTES, 0);
        this->m_dma_hls_fft_run_hdls[i].set_arg(0, this->m_dma_hls_fft_buffers[i]);
        this->m_dma_hls_fft_run_hdls[i].start();
        cols_buff_async_hdls.push_back(buffers_out[i].async("fftColsGraph[" + std::to_string(i) + "].gmio_out", 
                                                            XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                            BLOCK_SIZE_BYTES, 
                                                            offset));
    }
    
    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    for(int i = 0; i < num_of_buffers; i++) {
        cols_buff_async_hdls[i].wait();
        this->m_dma_hls_fft_run_hdls[i].wait();
    }
}

void ImgFFTCrossCorr::ifft2D(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
    int offset = 0;
    std::vector<xrt::bo::async_handle> cols_buff_async_hdls;
    std::vector<xrt::bo::async_handle> rows_buff_async_hdls;

    for(int i = 0; i < num_of_buffers; i++) {
        buffers_in[i].async("ifftColsGraph[" + std::to_string(i) + "].gmio_in", 
                            XCL_BO_SYNC_BO_GMIO_TO_AIE, BLOCK_SIZE_BYTES, offset);

        cols_buff_async_hdls.push_back(m_aie_to_pl_buffers[i].async("ifftColsGraph[" + std::to_string(i) + "].gmio_out", 
                                                                    XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                                    BLOCK_SIZE_BYTES, 
                                                                    offset));
    }

    // Setup and run PL kernels for data transfer to AIE
    for(int i = 0; i < num_of_buffers; i++) {
        cols_buff_async_hdls[i].wait();
        this->m_dma_hls_ifft_buffers[i] = xrt::bo(m_aie_to_pl_buffers[i]);
        this->m_dma_hls_ifft_buffers[i].sync(XCL_BO_SYNC_BO_TO_DEVICE, BLOCK_SIZE_BYTES, 0);
        this->m_dma_hls_ifft_run_hdls[i].set_arg(0, this->m_dma_hls_ifft_buffers[i]);
        this->m_dma_hls_ifft_run_hdls[i].start();
        rows_buff_async_hdls.push_back(buffers_out[i].async("ifftRowsGraph[" + std::to_string(i) + "].gmio_out", 
                                                            XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                            BLOCK_SIZE_BYTES, 
                                                            offset));
    }

    for(int i = 0; i < num_of_buffers; i++) {
        rows_buff_async_hdls[i].wait();
        this->m_dma_hls_ifft_run_hdls[i].wait();
    }
}

void ImgFFTCrossCorr::cplxConj(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
    int offset = 0;
    std::vector<xrt::bo::async_handle> buffer_async_hdls;

    for(int i = 0; i < num_of_buffers; i++) {
        buffers_in[i].async("cplxConjGraph[" + std::to_string(i) + "].gmio_in", 
                            XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                            BLOCK_SIZE_BYTES, 
                            offset);
        buffer_async_hdls.push_back(buffers_out[i].async("cplxConjGraph[" + std::to_string(i) + "].gmio_out",
                                                         XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                         BLOCK_SIZE_BYTES, 
                                                         offset));
    }

    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    for(int i = 0; i < num_of_buffers; i++) {
        buffer_async_hdls[i].wait();
    }
}

void ImgFFTCrossCorr::elemMatMult(xrt::aie::bo* buffersA_in, xrt::aie::bo* buffersB_in,
                                  xrt::aie::bo* buffers_out, int num_of_buffers) {
    std::vector<xrt::bo::async_handle> buffer_async_hdls;
    int per_ssr_byte_size = BLOCK_SIZE_BYTES / TP_SSR;

    for(int i = 0; i < num_of_buffers; i++) {
        for (int j = 0; j < TP_SSR; j++) {
            buffersA_in[i].async("hpGraph[" + std::to_string(i) + "].gmio_in_A[" + std::to_string(j) + "]", 
                                 XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                 per_ssr_byte_size, 
                                 j*per_ssr_byte_size);

            buffersB_in[i].async("hpGraph[" + std::to_string(i) + "].gmio_in_B[" + std::to_string(j) + "]", 
                                 XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                 per_ssr_byte_size, 
                                 j*per_ssr_byte_size);

            buffer_async_hdls.push_back(buffers_out[i].async("hpGraph[" + std::to_string(i) + "].gmio_out[" + std::to_string(j) + "]",
                                                             XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                             per_ssr_byte_size, 
                                                             j*per_ssr_byte_size));
        }
    }
    
    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    for(int i = 0; i < num_of_buffers; i++) {
        for (int j = 0; j < TP_SSR; j++) {
            buffer_async_hdls[i*TP_SSR + j].wait();
        }
    }
}

void ImgFFTCrossCorr::peakSearch(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
    int offset = 0;
    std::vector<xrt::bo::async_handle> buffer_async_hdls;

    for(int i = 0; i < num_of_buffers; i++) {
        buffers_in[i].async("peakGraph[" + std::to_string(i) + "].gmio_in", 
                            XCL_BO_SYNC_BO_GMIO_TO_AIE, BLOCK_SIZE_BYTES, offset);

        buffer_async_hdls.push_back(buffers_out[i].async("peakGraph[" + std::to_string(i) + "].gmio_out", 
                                                         XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                         sizeof(TT_DATA), 
                                                         offset));
    }

    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    // TODO: THIS wait() is not working as expected...maybe because I'm only asking for one byte?
    // Going to perform manual polling for data in calling code to see when it updated
    for(int i = 0; i < num_of_buffers; i++) {
        buffer_async_hdls[i].wait();
    }
}

ImgFFTCrossCorr::PeakValue ImgFFTCrossCorr::postPeakSearchOffset(TT_DATA peak_val) {
    PeakValue derived_peak;

    int idx = static_cast<int>(peak_val.imag);

    int row =  idx / MAT_COLS;
    int col = idx % MAT_ROWS;

    derived_peak.max.real = peak_val.real;
    derived_peak.x = (row + (MAT_ROWS/2)) % MAT_ROWS;
    derived_peak.y = (col + (MAT_COLS/2)) % MAT_COLS;

    return derived_peak;
}

void ImgFFTCrossCorr::printResults(PeakValue true_peak, PeakValue derived_peak) {
    const int x_diff = derived_peak.x - true_peak.x;
    const int y_diff = derived_peak.y - true_peak.y;
    const float euclid_dist = std::sqrt(std::pow(x_diff, 2) + std::pow(y_diff, 2));

    printf("\nPeak Search Info:\n");
    printf("peak.max[%d, %d] = %.1f\n", derived_peak.x, derived_peak.y, derived_peak.max.real);
    printf("x_error = %d\n", x_diff);
    printf("y_error = %d\n", y_diff);
    printf("euclid_dist = %.1f\n", euclid_dist);
    if(euclid_dist == 0) {
        printf("ITERATION PASSED!!!\n\n");
    } else {
        printf("ITERATION FAILED!!!\n\n");
    }
}
