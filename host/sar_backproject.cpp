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
, m_range_data_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_range_data_array(m_range_data_buffer.map<TT_DATA*>())
, m_ref_func_range_comp_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
, m_ref_func_range_comp_array(m_ref_func_range_comp_buffer.map<TT_DATA*>())
{
    //hid_t file_id;
    //file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //H5Fclose(file_id);

    // Instantiate kernels, handlers, and buffers for multiple instances
    for(int i = 0; i < this->m_instances; i++) {
        //std::string dma_hls_fft_kernel_str = "dma_hls:{dma_hls_" + std::to_string(i * 2) + "}";
        //std::string dma_hls_ifft_kernel_str = "dma_hls:{dma_hls_" + std::to_string(i * 2 + 1) + "}";
        std::string fftGraph_str = "fftGraph[" + std::to_string(i) + "]";
        std::string ifftGraph_str = "ifftGraph[" + std::to_string(i) + "]";
        //std::string cplxConjGraph_str = "cplxConjGraph[" + std::to_string(i) + "]";
        std::string hpGraph_str = "hpGraph[" + std::to_string(i) + "]";

        //m_dma_hls_fft_kernels.push_back(xrt::kernel(m_device, m_uuid, dma_hls_fft_kernel_str));
        //m_dma_hls_fft_run_hdls.push_back(xrt::run(m_dma_hls_fft_kernels[i]));
        //m_dma_hls_fft_buffers.push_back(xrt::bo(m_device, BLOCK_SIZE_BYTES, m_dma_hls_fft_kernels[i].group_id(0)));
        //m_dma_hls_ifft_kernels.push_back(xrt::kernel(m_device, m_uuid, dma_hls_ifft_kernel_str));
        //m_dma_hls_ifft_run_hdls.push_back(xrt::run(m_dma_hls_ifft_kernels[i]));
        //m_dma_hls_ifft_buffers.push_back(xrt::bo(m_device, BLOCK_SIZE_BYTES, m_dma_hls_ifft_kernels[i].group_id(0)));
        m_fft_graph_hdls.push_back(xrt::graph(m_device, m_uuid, fftGraph_str));
        m_ifft_graph_hdls.push_back(xrt::graph(m_device, m_uuid, ifftGraph_str));
        //m_cplx_conj_graph_hdls.push_back(xrt::graph(m_device, m_uuid, cplxConjGraph_str));
        m_hp_graph_hdls.push_back(xrt::graph(m_device, m_uuid, hpGraph_str));
        //m_aie_to_pl_buffers.push_back(xrt::aie::bo(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0));
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

void SARBackproject::runGraphs() {
    for(int i = 0; i < this->m_instances; i++) {
        this->m_fft_graph_hdls[i].run(MAT_ROWS*this->m_iter);
        this->m_ifft_graph_hdls[i].run(MAT_ROWS*this->m_iter);
        //this->m_cplx_conj_graph_hdls[i].run(MAT_ROWS*this->m_iter);
        this->m_hp_graph_hdls[i].run((MAT_ROWS/TP_NUM_FRAMES)*this->m_iter);
    }
}

void SARBackproject::endGraphs() {
    for(int i = 0; i < this->m_instances; i++) {
        this->m_fft_graph_hdls[i].end();
        this->m_ifft_graph_hdls[i].end();
        //this->m_cplx_conj_graph_hdls[i].end();
        this->m_hp_graph_hdls[i].end();
    }
}

bool SARBackproject::openHDF5File() {

    std::vector<TT_DATA> hh;
    std::vector<float> lut;

    hid_t file = H5Fopen(this->m_hdf5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (file < 0) {
        std::cerr << "Failed to open HDF5 file: " << this->m_hdf5_filename << std::endl;
        return false;
    }

    readCplxDataset(file, "science/LSAR/RRSD/swaths/frequencyA/txH/rxH/HH", hh);
    readFloatDataset(file, "science/LSAR/RRSD/swaths/frequencyA/txH/rxH/BFPQLUT", lut);
    
    for (int i=0; i<BLOCK_SIZE_ENTRIES; i++) {
        this->m_range_data_array[i] = (TT_DATA) {lut[hh[i].real], lut[hh[i].imag]};
    }

    for(int r=0; r<5; r++) {
        for(int c=0; c<5; c++) {
            int index = (r*MAT_COLS) + c;
            printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", r, c, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
        }
    }

    //TODO: Temporary until we know which ref function to use for range compression
    for(int r = 0; r < MAT_ROWS; r++) {
        for(int c = 0; c < MAT_COLS; c++) {
            int index = (r*MAT_COLS)+c;
            this->m_ref_func_range_comp_array[index] = (TT_DATA) {1, 0};
        }
    }


    return true;
}

void SARBackproject::fft(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {

    std::vector<xrt::bo::async_handle> buff_async_hdls;
    int ssr = (1 << TP_PARALLEL_POWER);
    int per_ssr_byte_size = BLOCK_SIZE_BYTES / ssr;

    //////////////////////TODO: TESTING...REMOVE WHEN DONE///////////////////////////// 
    //// Create a temporary array to hold the rearranged data
    //TT_DATA *temp = (TT_DATA *)malloc(BLOCK_SIZE_BYTES);
    //
    //// Rearrange the data into the temp array
    //int chunk_size = MAT_COLS;
    //int num_datasets = BLOCK_SIZE_ENTRIES / chunk_size;

    //for (int dataset = 0; dataset < num_datasets; dataset++) {
    //    int offset = dataset * chunk_size;

    //    for (int i = 0; i < chunk_size; i++) {
    //        int new_index = offset + (i / 4) + (i % 4) * (chunk_size / 4);
    //        temp[new_index] = this->m_range_data_array[offset + i];
    //    }
    //}

    //// Copy the rearranged data back into the original array
    //for (int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
    //    this->m_range_data_array[i] = temp[i];
    //}

    //for(int r=0; r<5; r++) {
    //    for(int c=0; c<5; c++) {
    //        int index = (r*MAT_COLS) + c;
    //        printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", r, c, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //    }
    //}
    //int index = (0*MAT_COLS) + 2047;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 2047, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 2048;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 2048, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 2049;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 2049, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 4095;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 4095, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 4096;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 4096, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 4097;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 4097, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);

    //// Free the temporary array
    //free(temp);
    //////////////////////TODO: TESTING...REMOVE WHEN DONE///////////////////////////// 

    for(int i = 0; i < num_of_buffers; i++) {
        for (int j = 0; j < ssr; j++) {
            buffers_in[i].async("fftGraph[" + std::to_string(i) + "].gmio_in[" + std::to_string(j) + "]", 
                                XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                per_ssr_byte_size, 
                                j*per_ssr_byte_size);

            buff_async_hdls.push_back(buffers_out[i].async("fftGraph[" + std::to_string(i) + "].gmio_out[" + std::to_string(j) + "]",
                                                           XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                           per_ssr_byte_size, 
                                                           j*per_ssr_byte_size));
        }
    }
    
    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    for(int i = 0; i < num_of_buffers; i++) {
        for (int j = 0; j < ssr; j++) {
            buff_async_hdls[i*ssr + j].wait();
        }
    }
}

void SARBackproject::ifft(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {

    std::vector<xrt::bo::async_handle> buff_async_hdls;
    int ssr = (1 << TP_PARALLEL_POWER);
    int per_ssr_byte_size = BLOCK_SIZE_BYTES / ssr;

    //////////////////////TODO: TESTING...REMOVE WHEN DONE///////////////////////////// 
    //// Create a temporary array to hold the rearranged data
    //TT_DATA *temp = (TT_DATA *)malloc(BLOCK_SIZE_BYTES);
    //
    //// Rearrange the data into the temp array
    //int chunk_size = MAT_COLS;
    //int num_datasets = BLOCK_SIZE_ENTRIES / chunk_size;

    //for (int dataset = 0; dataset < num_datasets; dataset++) {
    //    int offset = dataset * chunk_size;

    //    for (int i = 0; i < chunk_size; i++) {
    //        int new_index = offset + (i / 4) + (i % 4) * (chunk_size / 4);
    //        temp[new_index] = this->m_range_data_array[offset + i];
    //    }
    //}

    //// Copy the rearranged data back into the original array
    //for (int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
    //    this->m_range_data_array[i] = temp[i];
    //}

    //for(int r=0; r<5; r++) {
    //    for(int c=0; c<5; c++) {
    //        int index = (r*MAT_COLS) + c;
    //        printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", r, c, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //    }
    //}
    //int index = (0*MAT_COLS) + 2047;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 2047, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 2048;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 2048, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 2049;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 2049, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 4095;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 4095, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 4096;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 4096, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
    //index = (0*MAT_COLS) + 4097;
    //printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", 0, 4097, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);

    //// Free the temporary array
    //free(temp);
    //////////////////////TODO: TESTING...REMOVE WHEN DONE///////////////////////////// 
    

    for(int i = 0; i < num_of_buffers; i++) {
        for (int j = 0; j < ssr; j++) {
            buffers_in[i].async("ifftGraph[" + std::to_string(i) + "].gmio_in[" + std::to_string(j) + "]", 
                                XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                per_ssr_byte_size, 
                                j*per_ssr_byte_size);

            buff_async_hdls.push_back(buffers_out[i].async("ifftGraph[" + std::to_string(i) + "].gmio_out[" + std::to_string(j) + "]",
                                                           XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                                                           per_ssr_byte_size, 
                                                           j*per_ssr_byte_size));
        }
    }
    
    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    for(int i = 0; i < num_of_buffers; i++) {
        for (int j = 0; j < ssr; j++) {
            buff_async_hdls[i*ssr + j].wait();
        }
    }
}

void SARBackproject::cplxConj(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
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

void SARBackproject::elemMatMult(xrt::aie::bo* buffersA_in, xrt::aie::bo* buffersB_in,
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
