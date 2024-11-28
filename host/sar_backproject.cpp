// By: Austin Owens
// Date: 6/27/2024
// Desc: Class to handle SAR backprojection

#include "sar_backproject.h"
#include "../common.h"
#include <omp.h>
#include <fstream>
#include <thread>
#include <chrono>
#include <vector>
#include <utility>

double SARBackproject::total_time = 0;
double SARBackproject::total_avg_time = 0;
struct timespec SARBackproject::time_start = {0, 0};
struct timespec SARBackproject::time_end = {0, 0};

SARBackproject::SARBackproject(const char* xclbin_filename, const char* h5_filename, int iter, int instances)
: m_xclbin_filename(xclbin_filename)
, m_hdf5_filename(h5_filename)
, m_iter(iter)
, m_instances(instances)
, m_device(0)
, m_uuid(m_device.load_xclbin(this->m_xclbin_filename))
, m_x_ant_pos_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_x_ant_pos_array(m_x_ant_pos_buffer.map<float*>())
, m_y_ant_pos_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_y_ant_pos_array(m_y_ant_pos_buffer.map<float*>())
, m_z_ant_pos_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_z_ant_pos_array(m_z_ant_pos_buffer.map<float*>())
, m_ref_range_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_ref_range_array(m_ref_range_buffer.map<float*>())
, m_xy_px_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_xy_px_array(m_xy_px_buffer.map<TT_DATA*>())
, m_rc_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_rc_array(m_rc_buffer.map<TT_DATA*>())
, m_img_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_img_array(m_img_buffer.map<TT_DATA*>())
//, m_range_data_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
//, m_range_data_array(m_range_data_buffer.map<TT_DATA*>())
//, m_ref_func_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
//, m_ref_func_array(m_ref_func_buffer.map<TT_DATA*>())
{
    //hid_t file_id;
    //file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //H5Fclose(file_id);

    // Instantiate kernels, handlers, and buffers for multiple instances
    for(int i = 0; i < this->m_instances; i++) {
        std::string bpGraph_str = "bpGraph[" + std::to_string(i) + "]";
        ////std::string dma_hls_fft_kernel_str = "dma_hls:{dma_hls_" + std::to_string(i * 2) + "}";
        ////std::string dma_hls_ifft_kernel_str = "dma_hls:{dma_hls_" + std::to_string(i * 2 + 1) + "}";
        //std::string fftGraph_str = "fftGraph[" + std::to_string(i) + "]";
        //std::string ifftGraph_str = "ifftGraph[" + std::to_string(i) + "]";
        //std::string cplxConjGraph_str = "cplxConjGraph[" + std::to_string(i) + "]";
        //std::string hpGraph_str = "hpGraph[" + std::to_string(i) + "]";

        ////m_dma_hls_fft_kernels.push_back(xrt::kernel(m_device, m_uuid, dma_hls_fft_kernel_str));
        ////m_dma_hls_fft_run_hdls.push_back(xrt::run(m_dma_hls_fft_kernels[i]));
        ////m_dma_hls_fft_buffers.push_back(xrt::bo(m_device, BLOCK_SIZE_BYTES, m_dma_hls_fft_kernels[i].group_id(0)));
        ////m_dma_hls_ifft_kernels.push_back(xrt::kernel(m_device, m_uuid, dma_hls_ifft_kernel_str));
        ////m_dma_hls_ifft_run_hdls.push_back(xrt::run(m_dma_hls_ifft_kernels[i]));
        ////m_dma_hls_ifft_buffers.push_back(xrt::bo(m_device, BLOCK_SIZE_BYTES, m_dma_hls_ifft_kernels[i].group_id(0)));
        m_bp_graph_hdls.push_back(xrt::graph(m_device, m_uuid, bpGraph_str));
        //m_fft_graph_hdls.push_back(xrt::graph(m_device, m_uuid, fftGraph_str));
        //m_ifft_graph_hdls.push_back(xrt::graph(m_device, m_uuid, ifftGraph_str));
        //m_cplx_conj_graph_hdls.push_back(xrt::graph(m_device, m_uuid, cplxConjGraph_str));
        //m_hp_graph_hdls.push_back(xrt::graph(m_device, m_uuid, hpGraph_str));
        ////m_aie_to_pl_buffers.push_back(xrt::aie::bo(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0));
    }

    // Seed random number generator used for selecting
    //std::time_t seed = std::time(0);
    //std::time_t seed = 1719629552;
    //std::srand(seed);
    //std::cout << "Rand Seed: " << seed << std::endl;
}

void SARBackproject::startTime() {
    clock_gettime(CLOCK_MONOTONIC, &time_start);
}

void SARBackproject::endTime() {
    clock_gettime(CLOCK_MONOTONIC, &time_end);
}

void SARBackproject::resetTimer() {
    total_time = 0;
    total_avg_time = 0;
}

void SARBackproject::printTimeDiff(const char *msg) {
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

void SARBackproject::printTotalTime(int curr_iter) {
    printf("Elapsed total time: %.2f milliseconds\n", total_time);

    // Don't include the first iteration since it includes calculating the map image
    if(curr_iter > 0) {
        total_avg_time += total_time;
    }
    total_time = 0;
}

void SARBackproject::printAvgTime(int iterations) {
    if (iterations > 1) {
        total_avg_time = total_avg_time/(iterations-1);
        printf("Elapsed total average time over this last %d iteration(s): %.2f milliseconds\n", iterations-1, total_avg_time);
    }
}

void SARBackproject::reshapeMatrix(TT_DATA* input, int rows, int cols, int segment, bool reverse) {
    // Malloc tmp array 
    TT_DATA* tmp = (TT_DATA*) malloc(rows*cols*sizeof(TT_DATA));

    // Determine size of column segment
    int entry_seg = cols / segment;

    if (!reverse) {
        // Loop through each entry_seg
        #pragma omp parallel for collapse(2)
        for (int q = 0; q < segment; q++) {
            // Loop through each row
            for (int r = 0; r < rows; r++) {
                // Copy the corresponding part of each row into the tmp array
                memcpy(tmp + r * entry_seg + q * rows * entry_seg,
                       input + r * cols + q * entry_seg,
                       entry_seg * sizeof(TT_DATA));
            }
        }
    } else {
        #pragma omp parallel for collapse(2)
        for (int q = 0; q < 4; q++) {
            for (int r = 0; r < rows; r++) {
                memcpy(tmp + r * cols + q * entry_seg, 
                       input + r * entry_seg + q * rows * entry_seg, 
                       entry_seg * sizeof(TT_DATA));
            }
        }

    }

    // Copy the tmp array back into input
    memcpy(input, tmp, rows*cols*sizeof(TT_DATA));
    free(tmp);
}

void SARBackproject::strideCols(TT_DATA* input, int rows, int cols, int n, bool reverse) {
    // Allocate temporary array for rearranging
    TT_DATA* tmp = (TT_DATA*) malloc(rows*cols*sizeof(TT_DATA));
    int offset;
    
    // Rearrange (either forward or reverse)
    if (!reverse) {
        #pragma omp parallel for collapse(2)
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                offset = r * cols;
                int new_index = offset + (c / 4) + (c % 4) * (cols / 4);
                tmp[new_index] = input[offset + c];
            }
        }
    } else {
        #pragma omp parallel for collapse(2)
        for (int r = 0; r < rows; r++) {
            for (int c = 0; c < cols; c++) {
                offset = r * cols;
                int new_index = offset + (c / 4) + (c % 4) * (cols / 4);
                tmp[offset + c] = input[new_index];
            }
        }
    }

    // Copy rearranged data back into the original input array
    memcpy(input, tmp, rows*cols*sizeof(TT_DATA));

    // Free the temporary array
    free(tmp);
}

bool SARBackproject::readCplxDataset(hid_t file, const std::string& dataset_name, std::vector<TT_DATA>& data) {

    // Open the dataset
    hid_t dataset = H5Dopen(file, dataset_name.c_str(), H5P_DEFAULT);
    if (dataset < 0) {
        std::cerr << "Failed to open dataset: " << dataset_name << std::endl;
        return false;
    }

    // Get the dataspace and the dimensions of the dataset
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[2];
    H5Sget_simple_extent_dims(dataspace, dims, NULL); // Get dataset dimensions

    // Resize the 1D vector to hold the data
    data.resize(dims[0] * dims[1]);

    // Define the compound type for the TT_DATA structure
    hid_t compound_type = H5Tcreate(H5T_COMPOUND, sizeof(TT_DATA));
    H5Tinsert(compound_type, "r", HOFFSET(TT_DATA, real), H5T_NATIVE_FLOAT);
    H5Tinsert(compound_type, "i", HOFFSET(TT_DATA, imag), H5T_NATIVE_FLOAT);

    // Read the data into the vector
    H5Dread(dataset, compound_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

    // Close the resources
    H5Tclose(compound_type);
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return true;
}

bool SARBackproject::readFloatDataset(hid_t file, const std::string& dataset_name, std::vector<float>& data) {

    // Open the dataset
    hid_t dataset = H5Dopen(file, dataset_name.c_str(), H5P_DEFAULT);
    if (dataset < 0) {
        std::cerr << "Failed to open dataset: " << dataset_name << std::endl;
        return false;
    }

    // Get the dataspace and the dimensions of the dataset
    hid_t dataspace = H5Dget_space(dataset);
    hsize_t dims[1];
    H5Sget_simple_extent_dims(dataspace, dims, NULL);  // Get dataset dimensions

    // Resize the vector to hold the data
    data.resize(dims[0]);

    // Read the float32 data into the vector
    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());

    // Close the resources
    H5Sclose(dataspace);
    H5Dclose(dataset);

    return true;
}

bool SARBackproject::fetchRadarData() {

    std::vector<TT_DATA> hh;
    std::vector<float> lut;

    hid_t file = H5Fopen(this->m_hdf5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "Failed to open HDF5 file: " << this->m_hdf5_filename << std::endl;
        return false;
    }

    readCplxDataset(file, "science/LSAR/RRSD/swaths/frequencyA/txH/rxH/HH", hh);
    readFloatDataset(file, "science/LSAR/RRSD/swaths/frequencyA/txH/rxH/BFPQLUT", lut);
    
    //for (int i=0; i<BLOCK_SIZE_ENTRIES; i++) {
    //    this->m_range_data_array[i] = (TT_DATA) {lut[hh[i].real], lut[hh[i].imag]};
    //}

    const int PULSES = 1;
    for(int i = 0; i < PULSES; i++) {
        this->m_x_ant_pos_array[i] = 7089.2646;
        this->m_y_ant_pos_array[i] = 0.5289;
        this->m_z_ant_pos_array[i] = 7275.6719;
        this->m_ref_range_array[i] = 10158.399;
    }

    for(int i = 0; i < TP_POINT_SIZE; i++) {
        this->m_rc_array[i] = (TT_DATA) {i, i};
    }

    const float C = 299792458.0;
    float range_freq_step = 1471301.6;
    int half_range_samples = TP_POINT_SIZE/2;
    float min_freq = 9288080400.0;

    float range_width = C/(2.0*range_freq_step);
    float range_res = range_width/TP_POINT_SIZE;

    for(int i = 0; i < TP_POINT_SIZE; i++) {
        this->m_xy_px_array[i] = (TT_DATA) {(i-half_range_samples)*range_res, 0};
    }

    ////TODO: DEBUG
    //for(int r = 0; r < MAT_ROWS; r++) {
    //    for(int c = 0; c < MAT_COLS; c++) {
    //        int index = (r*MAT_COLS)+c;
    //        if (c == 1) {
    //            this->m_range_data_array[index] = (TT_DATA) {1.5, 1.5};
    //        }
    //        else {
    //            this->m_range_data_array[index] = (TT_DATA) {0, 0};
    //        }
    //    }
    //}

    //for(int r=0; r<5; r++) {
    //    for(int c=0; c<5; c++) {
    //        int index = (r*MAT_COLS) + c;
    //        printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", r, c, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //    }
    //}

    ////TODO: Temporary until we know which ref function to use for range compression
    //for(int r = 0; r < MAT_ROWS; r++) {
    //    for(int c = 0; c < MAT_COLS; c++) {
    //        int index = (r*MAT_COLS)+c;
    //        this->m_ref_func_array[index] = (TT_DATA) {1, 1};
    //    }
    //}

    return true;
}

void SARBackproject::runGraphs() {
    // Run graphs indefinitley (instead of for specific amount of iterations)
    for(int i = 0; i < this->m_instances; i++) {
        this->m_bp_graph_hdls[i].run(0);
        //this->m_fft_graph_hdls[i].run(0);
        //this->m_ifft_graph_hdls[i].run(0);
        //this->m_cplx_conj_graph_hdls[i].run(0);
        //this->m_hp_graph_hdls[i].run(0);
    }
}

void SARBackproject::bp(xrt::aie::bo* buffers_x_ant_pos_in, xrt::aie::bo* buffers_y_ant_pos_in, 
                        xrt::aie::bo* buffers_z_ant_pos_in, xrt::aie::bo* buffers_ref_range_in,
                        xrt::aie::bo* buffers_xy_px_in, xrt::aie::bo* buffers_rc_in, 
                        xrt::aie::bo* buffers_img_out, int num_of_buffers) {
    
    std::vector<xrt::bo::async_handle> buff_async_hdls;
    
    int pulses = 1;
    int per_bp_byte_size = (TP_POINT_SIZE/IMG_SOLVERS)*sizeof(TT_DATA);

    int32 rtp_valid_low_bound_result[4] = {};
    int32 rtp_valid_high_bound_result[4] = {};

    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
        for (int rc_seg = 0; rc_seg < IMG_SOLVERS; rc_seg++) {
            //printf("rc_seg = %d\n", rc_seg);
            
            // Data going to slowtime splicer
            buffers_x_ant_pos_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_x_ant_pos", 
                                                 XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                                 pulses*sizeof(float), 
                                                 0);
            buffers_y_ant_pos_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_y_ant_pos", 
                                                 XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                                 pulses*sizeof(float), 
                                                 0);
            buffers_z_ant_pos_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_z_ant_pos", 
                                                 XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                                 pulses*sizeof(float), 
                                                 0);
            buffers_ref_range_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_ref_range", 
                                                 XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                                 pulses*sizeof(float), 
                                                 0);

            for (int kern_id = 0; kern_id < IMG_SOLVERS; kern_id++) {
                //printf("kern_id = %d\n", kern_id);
                buffers_xy_px_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_xy_px[" + std::to_string(kern_id) + "]", 
                                                 XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                                 per_bp_byte_size, 
                                                 kern_id*per_bp_byte_size);

                buffers_rc_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_rc[" + std::to_string(kern_id) + "]", 
                                              XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                              per_bp_byte_size, 
                                              (3-rc_seg)*per_bp_byte_size);

                buff_async_hdls.push_back(buffers_img_out[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_out_img[" + std::to_string(kern_id) + "]",
                                                                          XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                                          per_bp_byte_size, 
                                                                          kern_id*per_bp_byte_size));
                
                this->m_bp_graph_hdls[buff_idx].read("bpGraph[" + std::to_string(buff_idx) + "].rtp_valid_low_bound[" + std::to_string(kern_id) + "]", 
                                                     rtp_valid_low_bound_result[kern_id]);
                //printf("Low Bound: %d\n", rtp_valid_low_bound_result[kern_id]);
                this->m_bp_graph_hdls[buff_idx].read("bpGraph[" + std::to_string(buff_idx) + "].rtp_valid_high_bound[" + std::to_string(kern_id) + "]", 
                                                     rtp_valid_high_bound_result[kern_id]);
                //printf("High Bound: %d\n", rtp_valid_high_bound_result[kern_id]);
            }


            for(int i=0; i<4; i++) {
                buff_async_hdls[i].wait();
                //for(int j=i*2048; j<(i*2048)+8; j++) {
                //    printf("%d: %f\t%f\n", j, this->m_img_array[j].real, this->m_img_array[j].imag);
                //}
            }
            buff_async_hdls.clear();
        }
    }

    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    //for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
    //    printf("Waiting on idx%d\n", idx);
    //    buff_async_hdls[idx].wait();
    //}
    //std::this_thread::sleep_for(std::chrono::seconds(1));
}

void SARBackproject::fft(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
    
    /* 
     * The FFT operation uses 2^TP_PARALLEL_POWER subframe processors in parallel each performing 
     * an FFT of N/(2^TP_PARALLEL_POWER) point size. These subframe outputs are then combined by 
     * TP_PARALLEL_POWER stages of radix2 to create the final result. As a result, the FFT input
     * requires the row of data to be split into each port. 
     *
     * For example, for 8192 point size and a TP_PARALLEL_POWER = 2, the size for each port
     * becomes 8192/4=2048. This means: 
     *
     * N[0-2047]    -> port0
     * N[2048-4095] -> port1
     * N[4096-6143] -> port2 
     * N[6144-8191] -> port3
     *
     * Because of the radix stages, the output is also mixed. Continuing with the example
     * above, assume the FFT output provides an array of P[8192]. To get the actual order
     * of the FFT (Q), the array would need to be rearranged as follows:
     *
     * P[0]  -> Q[0]
     * P[4]  -> Q[1]
     * P[8]  -> Q[2]
     * P[12] -> Q[3]
     * .
     * .
     * .
     * P[1]  -> Q[2048]
     * P[5]  -> Q[2049]
     * P[9]  -> Q[2050]
     * P[13] -> Q[2051]
     * .
     * .
     * .
     * P[2]  -> Q[4096]
     * P[6]  -> Q[4097]
     * P[10] -> Q[4098]
     * P[14] -> Q[4099]
     * .
     * .
     * .
     * P[3]  -> Q[6144]
     * P[7]  -> Q[6145]
     * P[11] -> Q[6146]
     * P[15] -> Q[6147]
     *
     * The reshapeMatrix() function can be used for reshaping the data for the input ports
     * and strideCols() function can be used for reshaping the data for the output ports.
     */

    std::vector<xrt::bo::async_handle> buff_async_hdls;

    // 54 MS TO COMPLETE
    int per_ssr_byte_size = BLOCK_SIZE_BYTES / FFT_NPORTS;

    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
        for (int ssr = 0; ssr < FFT_NPORTS; ssr++) {
            buffers_in[buff_idx].async("fftGraph[" + std::to_string(buff_idx) + "].gmio_in[" + std::to_string(ssr) + "]", 
                                       XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                       per_ssr_byte_size, 
                                       ssr*per_ssr_byte_size);

            buff_async_hdls.push_back(buffers_out[buff_idx].async("fftGraph[" + std::to_string(buff_idx) + "].gmio_out[" + std::to_string(ssr) + "]",
                                                                  XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                                  per_ssr_byte_size, 
                                                                  ssr*per_ssr_byte_size));
        }
    }

    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
        buff_async_hdls[idx].wait();
    }
}

void SARBackproject::ifft(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {

    std::vector<xrt::bo::async_handle> buff_async_hdls;
    int per_ssr_byte_size = BLOCK_SIZE_BYTES / FFT_NPORTS;

    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
        for (int ssr = 0; ssr < FFT_NPORTS; ssr++) {
            buffers_in[buff_idx].async("ifftGraph[" + std::to_string(buff_idx) + "].gmio_in[" + std::to_string(ssr) + "]", 
                                       XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                       per_ssr_byte_size, 
                                       ssr*per_ssr_byte_size);

            buff_async_hdls.push_back(buffers_out[buff_idx].async("ifftGraph[" + std::to_string(buff_idx) + "].gmio_out[" + std::to_string(ssr) + "]",
                                                                  XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                                  per_ssr_byte_size, 
                                                                  ssr*per_ssr_byte_size));
        }
    }

    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
        buff_async_hdls[idx].wait();
    }
}

void SARBackproject::cplxConj(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
    int offset = 0;
    std::vector<xrt::bo::async_handle> buff_async_hdls;

    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
        buffers_in[buff_idx].async("cplxConjGraph[" + std::to_string(buff_idx) + "].gmio_in", 
                                   XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                   BLOCK_SIZE_BYTES, 
                                   offset);
        buff_async_hdls.push_back(buffers_out[buff_idx].async("cplxConjGraph[" + std::to_string(buff_idx) + "].gmio_out",
                                                              XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                              BLOCK_SIZE_BYTES, 
                                                              offset));
    }

    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
        buff_async_hdls[idx].wait();
    }
}

void SARBackproject::elemMatMult(xrt::aie::bo* buffersA_in, xrt::aie::bo* buffersB_in,
                                  xrt::aie::bo* buffers_out, int num_of_buffers) {
    std::vector<xrt::bo::async_handle> buff_async_hdls;
    int per_ssr_byte_size = BLOCK_SIZE_BYTES / TP_SSR;

    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
        for (int ssr = 0; ssr < TP_SSR; ssr++) {
            buffersA_in[buff_idx].async("hpGraph[" + std::to_string(buff_idx) + "].gmio_in_A[" + std::to_string(ssr) + "]", 
                                        XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                        per_ssr_byte_size, 
                                        ssr*per_ssr_byte_size);

            buffersB_in[buff_idx].async("hpGraph[" + std::to_string(buff_idx) + "].gmio_in_B[" + std::to_string(ssr) + "]", 
                                        XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                        per_ssr_byte_size, 
                                        ssr*per_ssr_byte_size);

            buff_async_hdls.push_back(buffers_out[buff_idx].async("hpGraph[" + std::to_string(buff_idx) + "].gmio_out[" + std::to_string(ssr) + "]",
                                                                  XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                                  per_ssr_byte_size, 
                                                                  ssr*per_ssr_byte_size));
        }
    }
    
    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
        buff_async_hdls[idx].wait();
    }
}
