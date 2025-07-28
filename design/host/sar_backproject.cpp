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
#include <regex>

double SARBackproject::total_time = 0;
double SARBackproject::total_avg_time = 0;
struct timespec SARBackproject::time_start = {0, 0};
struct timespec SARBackproject::time_end = {0, 0};

SARBackproject::SARBackproject(const char* xclbin_filename, 
                               const char* st_dataset_filename, 
                               const char* rc_dataset_filename, 
                               const char* img_out_filename, 
                               int iter, 
                               int instances)
: m_xclbin_filename(xclbin_filename)
//, m_hdf5_filename(h5_filename)
, m_st_dataset_filename(st_dataset_filename)
, m_rc_dataset_filename(rc_dataset_filename)
, m_img_out_filename(img_out_filename)
, m_iter(iter)
, m_instances(instances)
, m_device(0)
, m_uuid(m_device.load_xclbin(this->m_xclbin_filename))
, m_broadcast_data_buffer(m_device, PULSES*BC_ELEMENTS*sizeof(float), xrt::bo::flags::normal, 0)
, m_broadcast_data_array(m_broadcast_data_buffer.map<float*>())
//, m_x_ant_pos_buffer(m_device, PULSES*sizeof(float), xrt::bo::flags::normal, 0)
//, m_x_ant_pos_array(m_x_ant_pos_buffer.map<float*>())
//, m_y_ant_pos_buffer(m_device, PULSES*sizeof(float), xrt::bo::flags::normal, 0)
//, m_y_ant_pos_array(m_y_ant_pos_buffer.map<float*>())
//, m_z_ant_pos_buffer(m_device, PULSES*sizeof(float), xrt::bo::flags::normal, 0)
//, m_z_ant_pos_array(m_z_ant_pos_buffer.map<float*>())
//, m_ref_range_buffer(m_device, PULSES*sizeof(float), xrt::bo::flags::normal, 0)
//, m_ref_range_array(m_ref_range_buffer.map<float*>())
, m_xyz_px_buffer(m_device, PULSES*RC_SAMPLES*sizeof(float)*3, xrt::bo::flags::normal, 0)
, m_xyz_px_array(m_xyz_px_buffer.map<float*>())
, m_rc_buffer(m_device, PULSES*RC_SAMPLES*sizeof(cfloat), xrt::bo::flags::normal, 0)
, m_rc_array(m_rc_buffer.map<cfloat*>())
//, m_img_buffer(m_device, PULSES*RC_SAMPLES*sizeof(cfloat), xrt::bo::flags::normal, 0)
//, m_img_array(m_img_buffer.map<cfloat*>())
//, m_range_data_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
//, m_range_data_array(m_range_data_buffer.map<cfloat*>())
//, m_ref_func_buffer(m_device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, 0)
//, m_ref_func_array(m_ref_func_buffer.map<cfloat*>())
{
    //hid_t file_id;
    //file_id = H5Fcreate("test.h5", H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);
    //H5Fclose(file_id);


    // Instantiate kernels, handlers, and buffers for multiple instances
    // TODO: Instances wont work right if it is greater than 1. This is fine
    // for now since there are not plans to instantiate the design more than
    // once, but would be nice if this still worked for future proofing
    for(int i = 0; i < this->m_instances; i++) {
        for(int sw_id = 0; sw_id < AIE_SWITCHES; sw_id++) {
            std::string dma_pkt_router_kernel_str = "dma_pkt_router:{dma_pkt_router_" + std::to_string(sw_id) + "}";
            m_dma_pkt_router_kernels.push_back(xrt::kernel(m_device, m_uuid, dma_pkt_router_kernel_str));
            m_dma_pkt_router_run_hdls.push_back(xrt::run(m_dma_pkt_router_kernels[sw_id]));
        }
        
        // All kernel instances map to the same DDR memory bank. This means
        // kernel1.group_id(1) should return the same bank ID as
        // kernel0.group_id(1). Therefore, it doesnt matter which instances
        // group_id is used, the result is the same memory bank.
        m_img_buffers.push_back(xrt::bo(m_device, PULSES*RC_SAMPLES*8, m_dma_pkt_router_kernels[0].group_id(1)));
        m_img_arrays.push_back(m_img_buffers[i].map<cfloat*>());

        std::string bpGraph_str = "bpGraph[" + std::to_string(i) + "]";
        m_bp_graph_hdls.push_back(xrt::graph(this->m_device, this->m_uuid, bpGraph_str));
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

void SARBackproject::unwrap(double* angles) {
    // Store the first original angle for difference calculations.
    double prev_orig = angles[0];
    // The first angle remains unchanged.
    for (int i = 1; i < PULSES; i++) {
        // Save the current original value before modifying it.
        double current_orig = angles[i];
        double diff = current_orig - prev_orig;
        // Wrap diff into the range [-pi, pi)
        double dp = fmod(diff + PI, TWO_PI);
        if (dp < 0)
            dp += TWO_PI;
        dp -= PI;
        // Adjust the edge case: if dp == -pi and diff > 0, set dp to pi.
        if (dp == -PI && diff > 0)
            dp = PI;
        // The new (unwrapped) angle is the previous unwrapped angle plus dp.
        angles[i] = angles[i - 1] + dp;
        // Update prev_orig to the original value before modification.
        prev_orig = current_orig;
    }
}

//void SARBackproject::reshapeMatrix(cfloat* input, int rows, int cols, int segment, bool reverse) {
//    // Malloc tmp array 
//    cfloat* tmp = (cfloat*) malloc(rows*cols*sizeof(cfloat));
//
//    // Determine size of column segment
//    int entry_seg = cols / segment;
//
//    if (!reverse) {
//        // Loop through each entry_seg
//        #pragma omp parallel for collapse(2)
//        for (int q = 0; q < segment; q++) {
//            // Loop through each row
//            for (int r = 0; r < rows; r++) {
//                // Copy the corresponding part of each row into the tmp array
//                memcpy(tmp + r * entry_seg + q * rows * entry_seg,
//                       input + r * cols + q * entry_seg,
//                       entry_seg * sizeof(cfloat));
//            }
//        }
//    } else {
//        #pragma omp parallel for collapse(2)
//        for (int q = 0; q < 4; q++) {
//            for (int r = 0; r < rows; r++) {
//                memcpy(tmp + r * cols + q * entry_seg, 
//                       input + r * entry_seg + q * rows * entry_seg, 
//                       entry_seg * sizeof(cfloat));
//            }
//        }
//
//    }
//
//    // Copy the tmp array back into input
//    memcpy(input, tmp, rows*cols*sizeof(cfloat));
//    free(tmp);
//}

//void SARBackproject::strideCols(cfloat* input, int rows, int cols, int n, bool reverse) {
//    // Allocate temporary array for rearranging
//    cfloat* tmp = (cfloat*) malloc(rows*cols*sizeof(cfloat));
//    int offset;
//    
//    // Rearrange (either forward or reverse)
//    if (!reverse) {
//        #pragma omp parallel for collapse(2)
//        for (int r = 0; r < rows; r++) {
//            for (int c = 0; c < cols; c++) {
//                offset = r * cols;
//                int new_index = offset + (c / 4) + (c % 4) * (cols / 4);
//                tmp[new_index] = input[offset + c];
//            }
//        }
//    } else {
//        #pragma omp parallel for collapse(2)
//        for (int r = 0; r < rows; r++) {
//            for (int c = 0; c < cols; c++) {
//                offset = r * cols;
//                int new_index = offset + (c / 4) + (c % 4) * (cols / 4);
//                tmp[offset + c] = input[new_index];
//            }
//        }
//    }
//
//    // Copy rearranged data back into the original input array
//    memcpy(input, tmp, rows*cols*sizeof(cfloat));
//
//    // Free the temporary array
//    free(tmp);
//}

//bool SARBackproject::readCplxDataset(hid_t file, const std::string& dataset_name, std::vector<cfloat>& data) {
//
//    // Open the dataset
//    hid_t dataset = H5Dopen(file, dataset_name.c_str(), H5P_DEFAULT);
//    if (dataset < 0) {
//        std::cerr << "Failed to open dataset: " << dataset_name << std::endl;
//        return false;
//    }
//
//    // Get the dataspace and the dimensions of the dataset
//    hid_t dataspace = H5Dget_space(dataset);
//    hsize_t dims[2];
//    H5Sget_simple_extent_dims(dataspace, dims, NULL); // Get dataset dimensions
//
//    // Resize the 1D vector to hold the data
//    data.resize(dims[0] * dims[1]);
//
//    // Define the compound type for the cfloat structure
//    hid_t compound_type = H5Tcreate(H5T_COMPOUND, sizeof(cfloat));
//    H5Tinsert(compound_type, "r", HOFFSET(cfloat, real), H5T_NATIVE_FLOAT);
//    H5Tinsert(compound_type, "i", HOFFSET(cfloat, imag), H5T_NATIVE_FLOAT);
//
//    // Read the data into the vector
//    H5Dread(dataset, compound_type, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
//
//    // Close the resources
//    H5Tclose(compound_type);
//    H5Sclose(dataspace);
//    H5Dclose(dataset);
//
//    return true;
//}
//
//bool SARBackproject::readFloatDataset(hid_t file, const std::string& dataset_name, std::vector<float>& data) {
//
//    // Open the dataset
//    hid_t dataset = H5Dopen(file, dataset_name.c_str(), H5P_DEFAULT);
//    if (dataset < 0) {
//        std::cerr << "Failed to open dataset: " << dataset_name << std::endl;
//        return false;
//    }
//
//    // Get the dataspace and the dimensions of the dataset
//    hid_t dataspace = H5Dget_space(dataset);
//    hsize_t dims[1];
//    H5Sget_simple_extent_dims(dataspace, dims, NULL);  // Get dataset dimensions
//
//    // Resize the vector to hold the data
//    data.resize(dims[0]);
//
//    // Read the float32 data into the vector
//    H5Dread(dataset, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, data.data());
//
//    // Close the resources
//    H5Sclose(dataspace);
//    H5Dclose(dataset);
//
//    return true;
//}
//
//bool SARBackproject::fetchRadarData() {
//
//    std::vector<cfloat> hh;
//    std::vector<float> lut;
//
//    hid_t file = H5Fopen(this->m_hdf5_filename, H5F_ACC_RDONLY, H5P_DEFAULT);
//    if (file < 0) {
//        std::cerr << "Failed to open HDF5 file: " << this->m_hdf5_filename << std::endl;
//        return false;
//    }
//
//    readCplxDataset(file, "science/LSAR/RRSD/swaths/frequencyA/txH/rxH/HH", hh);
//    readFloatDataset(file, "science/LSAR/RRSD/swaths/frequencyA/txH/rxH/BFPQLUT", lut);
//    
//    //for (int i=0; i<BLOCK_SIZE_ENTRIES; i++) {
//    //    this->m_range_data_array[i] = (cfloat) {lut[hh[i].real], lut[hh[i].imag]};
//    //}
//
//    const int PULSES = 1;
//    for(int i = 0; i < PULSES; i++) {
//        this->m_x_ant_pos_array[i] = 7089.2646;
//        this->m_y_ant_pos_array[i] = 0.5289;
//        this->m_z_ant_pos_array[i] = 7275.6719;
//        this->m_ref_range_array[i] = 10158.399;
//    }
//
//    for(int i = 0; i < RC_SAMPLES; i++) {
//        this->m_rc_array[i] = (cfloat) {i, i};
//    }
//
//    const float C = 299792458.0;
//    float range_freq_step = 1471301.6;
//    int half_range_samples = RC_SAMPLES/2;
//    float min_freq = 9288080400.0;
//
//    float range_width = C/(2.0*range_freq_step);
//    float range_res = range_width/RC_SAMPLES;
//
//    //for(int i = 0; i < RC_SAMPLES; i++) {
//    //    this->m_xy_px_array[i] = (cfloat) {(i-half_range_samples)*range_res, 0};
//    //}
//
//    ////TODO: DEBUG
//    //for(int r = 0; r < MAT_ROWS; r++) {
//    //    for(int c = 0; c < MAT_COLS; c++) {
//    //        int index = (r*MAT_COLS)+c;
//    //        if (c == 1) {
//    //            this->m_range_data_array[index] = (cfloat) {1.5, 1.5};
//    //        }
//    //        else {
//    //            this->m_range_data_array[index] = (cfloat) {0, 0};
//    //        }
//    //    }
//    //}
//
//    //for(int r=0; r<5; r++) {
//    //    for(int c=0; c<5; c++) {
//    //        int index = (r*MAT_COLS) + c;
//    //        printf("this->m_range_data_array[%d][%d] = {%f, %f}\n", r, c, this->m_range_data_array[index].real, this->m_range_data_array[index].imag);
//    //    }
//    //}
//
//    ////TODO: Temporary until we know which ref function to use for range compression
//    //for(int r = 0; r < MAT_ROWS; r++) {
//    //    for(int c = 0; c < MAT_COLS; c++) {
//    //        int index = (r*MAT_COLS)+c;
//    //        this->m_ref_func_array[index] = (cfloat) {1, 1};
//    //    }
//    //}
//
//    return true;
//}

bool SARBackproject::writeImg() {
    FILE *img_fp = fopen(this->m_img_out_filename, "w");
    if (img_fp == NULL) {
        std::cerr << "Error opening image output file!" << std::endl;
        return 1;
    }
    
    // Only looking at the first instance of m_img_arrays (assumes the design
    // is only instantiated once)
    fprintf(img_fp, "%.12f%+.12fi", this->m_img_arrays[0][0].real, this->m_img_arrays[0][0].imag);
    for(int i=1; i<PULSES*RC_SAMPLES; i++) {
        if (i%RC_SAMPLES == 0) {
            fprintf(img_fp, "\n");
        }
        fprintf(img_fp, ",%.12f%+.12fi", this->m_img_arrays[0][i].real, this->m_img_arrays[0][i].imag);
    }
    return 0;
}

bool SARBackproject::fetchRadarData() {

    // SLOWTIME DATA
    std::ifstream st_file(this->m_st_dataset_filename);
    if (!st_file.is_open()) {
        std::cerr << "Error opening slowtime dataset!" << std::endl;
        return 1;
    }
    std::string line;
    int pulse_idx = 0;
    while (std::getline(st_file, line) && pulse_idx < PULSES) {
        std::stringstream ss(line);
        std::string value;

        // Read each column from the line and store in respective array
        std::getline(ss, value, ',');
        this->m_broadcast_data_array[BC_ELEMENTS*pulse_idx] = std::stof(value);

        std::getline(ss, value, ',');
        this->m_broadcast_data_array[BC_ELEMENTS*pulse_idx + 1] = std::stof(value);

        std::getline(ss, value, ',');
        this->m_broadcast_data_array[BC_ELEMENTS*pulse_idx + 2] = std::stof(value);

        std::getline(ss, value, ',');
        this->m_broadcast_data_array[BC_ELEMENTS*pulse_idx + 3] = std::stof(value);

        pulse_idx++;
    }
    
    //for(int i=0; i<2; i++) {
    //    printf("this->m_x_ant_pos_array[%d] = %f\n", i, this->m_x_ant_pos_array[i]);
    //    printf("this->m_y_ant_pos_array[%d] = %f\n", i, this->m_y_ant_pos_array[i]);
    //    printf("this->m_z_ant_pos_array[%d] = %f\n", i, this->m_z_ant_pos_array[i]);
    //    printf("this->m_ref_range_array[%d] = %f\n", i, this->m_ref_range_array[i]);
    //}
    
    // RANGE COMPRESSED DATA
    std::ifstream rc_file(this->m_rc_dataset_filename);
    if (!rc_file.is_open()) {
        std::cerr << "Error opening range compressed dataset!" << std::endl;
        return 1;
    }
    pulse_idx = 0;
    while (std::getline(rc_file, line) && pulse_idx < PULSES) {
        std::stringstream ss(line);
        std::string value;
        int rc_samp_cnt = 0;

        // Read each column (complex number) from the line and store in array
        while (std::getline(ss, value, ',') && rc_samp_cnt < RC_SAMPLES) {
            std::regex complex_regex(R"(([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)([+-](?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)i)");
            std::smatch match;
            std::regex_search(value, match, complex_regex);

            float real_part = std::stof(match[1].str());
            float imag_part = std::stof(match[2].str());

            // Store the complex number in the array
            this->m_rc_array[pulse_idx*RC_SAMPLES + rc_samp_cnt] = (cfloat) {real_part, imag_part};

            rc_samp_cnt++;
        }
        rc_samp_cnt = 0;
        pulse_idx++;
    }

    return 0;
}

void SARBackproject::genTargetPixels() {
    // AZIMUTH RESOLUTION GRID
    double az_res = 0;
    double half_az_width = 0;
    if (PULSES != 1) {
        double az_ant[PULSES];
        for (int i = 0; i < PULSES; i++) {
            // atan2(y_ant_pos / x_ant_pos)
            az_ant[i] = std::atan2(this->m_broadcast_data_array[BC_ELEMENTS*i + 1], this->m_broadcast_data_array[BC_ELEMENTS*i]);
        }
        this->unwrap(az_ant);
        double sum_diff = 0.0;
        for (int i = 1; i < PULSES; i++) {
            sum_diff += (az_ant[i] - az_ant[i - 1]);
        }
        double mean_diff = sum_diff / (PULSES - 1);
        double delta_az = std::fabs(mean_diff);
        double min_az = *std::min_element(az_ant, az_ant + PULSES);
        double max_az = *std::max_element(az_ant, az_ant + PULSES);
        double total_az = max_az - min_az;
        az_res = C/(2.0*total_az*MIN_FREQ);
        double az_width = C/(2.0*delta_az*MIN_FREQ);
        half_az_width = az_width/2.0;
    }

    // TARGET PIXELS
    //for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
    //    for(int rng_idx = 0; rng_idx < RC_SAMPLES; rng_idx++) {
    //        int idx = pulse_idx*RC_SAMPLES + rng_idx;
    //        this->m_xy_px_array[idx] = (cfloat) {(rng_idx-HALF_RANGE_SAMPLES)*RANGE_RES, az_res*pulse_idx - half_az_width};
    //        this->m_z_px_array[idx] = 0.0;

    //        //printf("pixels[%d] = {%f, %f}\n", idx, this->m_xy_px_array[idx].real, this->m_xy_px_array[idx].imag);
    //    }
    //}
    int idx = 0;
    for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
        for(int rng_idx = 0; rng_idx < RC_SAMPLES; rng_idx++) {

            // X target pixels
            this->m_xyz_px_array[idx++] = (rng_idx-HALF_RANGE_SAMPLES)*RANGE_RES;

            // Y target pixels
            this->m_xyz_px_array[idx++] = az_res*pulse_idx - half_az_width;

            // Z target pixels
            this->m_xyz_px_array[idx++] = 0.0;

            //printf("pixels[%d] = {%f, %f, %f}\n", (idx-3)/3, this->m_xyz_px_array[idx-3], this->m_xyz_px_array[idx-2], this->m_xyz_px_array[idx-1]);
        }
    }
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

void SARBackproject::bp(xrt::aie::bo* buffers_broadcast_data_in, xrt::aie::bo* buffers_rc_in, 
        xrt::aie::bo* buffers_xyz_px_in, xrt::bo* buffers_img_out, int num_of_buffers) {
    
    std::vector<xrt::bo::async_handle> buff_async_hdls;
    
    // Number of target pixels for each px_demux_kern to process
    int px_per_demux_kern = (PULSES*RC_SAMPLES)/AIE_SWITCHES;
    //int px_per_ai = (PULSES*RC_SAMPLES)/IMG_SOLVERS;

    //int32 rtp_valid_low_bound_result[4] = {};
    //int32 rtp_valid_high_bound_result[4] = {};

    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
        // Slowtime and range compressed data going to data broadcaster splicer
        buffers_broadcast_data_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_st", 
                                                  XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                                  PULSES*BC_ELEMENTS*sizeof(float), 
                                                  0);
        for (int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
            buffers_rc_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_rc", 
                                          XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                          RC_SAMPLES*sizeof(cfloat), 
                                          (pulse_idx*RC_SAMPLES)*sizeof(cfloat));
            

            for (int sw_id=0; sw_id<AIE_SWITCHES; sw_id++) {
                buffers_xyz_px_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].bpCluster[" + std::to_string(sw_id) + "].gmio_in_xyz_px", 
                                                  XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                                                  px_per_demux_kern*sizeof(float)*3, 
                                                  sw_id*px_per_demux_kern*sizeof(float)*3);
            }

            for (int kern_id = 0; kern_id < IMG_SOLVERS; kern_id++) {
                // Dump image if on last pulse, otherwise keep focusing the image
                if (pulse_idx == PULSES-1) {
                    this->m_bp_graph_hdls[buff_idx].update("bpGraph[" + std::to_string(buff_idx) + "].rtp_dump_img_in[" + std::to_string(kern_id) + "]", 1);
                } else {
                    this->m_bp_graph_hdls[buff_idx].update("bpGraph[" + std::to_string(buff_idx) + "].rtp_dump_img_in[" + std::to_string(kern_id) + "]", 0);
                }

                //printf("kern_id = %d\n", kern_id);
                //buffers_xy_px_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_xy_px[" + std::to_string(kern_id) + "]", 
                //                                 XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                //                                 px_per_ai*sizeof(cfloat), 
                //                                 kern_id*px_per_ai*sizeof(cfloat));

                //buffers_z_px_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_z_px[" + std::to_string(kern_id) + "]", 
                //                                 XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                //                                 px_per_ai*sizeof(float), 
                //                                 kern_id*px_per_ai*sizeof(float));


                // Derive RC offsets for AI kernel per pulse. 
                // The motivation behind the following derivation is to try and pass "usable" range 
                // compressed values into each img_reconstruct_kern AI kernel. Because we calculate the 
                // boundaries of what the AI kernels will use here (on host/ARM), that gives us insight
                // into which range compressed values will actually be utilized for a specific AI kernel 
                // given their target pixels they are deriving for. This reduces the necessary size of the
                // rc_in ping-pong buffer into the kernel (which is necessary because we are already pushing
                // the stack limit). Note: This same calculation occurs on the AI kernel for each target pixel.
                //float dR_bounds = sqrt(pow(this->m_x_ant_pos_array[pulse_idx] - this->m_xy_px_array[kern_id*px_per_ai + (px_per_ai-1)].real, 2) 
                //        + pow(this->m_y_ant_pos_array[pulse_idx] - this->m_xy_px_array[kern_id*px_per_ai + (px_per_ai-1)].imag, 2) 
                //        + pow(this->m_z_ant_pos_array[pulse_idx] - this->m_z_px_array[kern_id*px_per_ai + (px_per_ai-1)], 2))
                //        - this->m_ref_range_array[pulse_idx];

                //float px_idx_bound = dR_bounds/RANGE_RES + HALF_RANGE_SAMPLES;
                //int rounded_px_idx_bound = (int) std::floor(px_idx_bound);

                //// The offset for gm2aie needs to be 128 bit aligned. Because cfloats are 8B (64b),
                //// then the rc_idx_offset just needs to be even
                //int rc_idx_offset = rounded_px_idx_bound - rounded_px_idx_bound%2;
                
                // TODO: Probably remove this
                //this->m_bp_graph_hdls[buff_idx].update("bpGraph[" + std::to_string(buff_idx) + "].rtp_rc_idx_offset_in[" + std::to_string(kern_id) + "]", 0);


                // Need to be wiser and stratigically pass in data based on what the input target pixels are for that AI tile
                //bpGraph[inst].gmio_in_rc[kern_id].gm2aie_nb(rc_array, RC_SAMPLES*sizeof(cfloat));
                //bpGraph[inst].gmio_in_rc[kern_id].gm2aie_nb(rc_array + (3-kern_id)*px_per_ai, px_per_ai*sizeof(cfloat));

                //buffers_rc_in[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_in_rc[" + std::to_string(kern_id) + "]", 
                //                              XCL_BO_SYNC_BO_GMIO_TO_AIE, 
                //                              RC_SAMPLES*sizeof(cfloat), 
                //                              (pulse_idx*RC_SAMPLES)*sizeof(cfloat));

                
                //buff_async_hdls.push_back(buffers_img_out[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_out_img[" + std::to_string(kern_id) + "]",
                //                                                          XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                //                                                          px_per_ai*sizeof(cfloat), 
                //                                                          kern_id*px_per_ai*sizeof(cfloat)));
                //buffers_img_out[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_out_img[" + std::to_string(kern_id) + "]",
                //                                XCL_BO_SYNC_BO_AIE_TO_GMIO, 
                //                                px_per_ai*sizeof(cfloat), 
                //                                kern_id*px_per_ai*sizeof(cfloat));

            }
            //buffers_img_out[buff_idx].async("bpGraph[" + std::to_string(buff_idx) + "].gmio_out_img[" + std::to_string(0) + "]",
            //                                XCL_BO_SYNC_BO_AIE_TO_GMIO, 
            //                                PULSES*RC_SAMPLES*sizeof(cfloat), 
            //                                0);




            //cols_buff_async_hdls.push_back(buffers_out[i].async("fftColsGraph[" + std::to_string(i) + "].gmio_out", 
            //                                                    XCL_BO_SYNC_BO_AIE_TO_GMIO, 
            //                                                    PULSES*RC_SAMPLES*sizeof(cfloat), 
            //                                                    offset));


            // Setup and run PL kernels for data transfer to AIE
            //for(int i = 0; i < num_of_buffers; i++) {
            //    rows_buff_async_hdls[i].wait();
            //    this->m_dma_hls_fft_buffers[i] = xrt::bo(m_aie_to_pl_buffers[i]);
            //    this->m_dma_hls_fft_buffers[i].sync(XCL_BO_SYNC_BO_TO_DEVICE, BLOCK_SIZE_BYTES, 0);
            //    this->m_dma_hls_fft_run_hdls[i].set_arg(0, this->m_dma_hls_fft_buffers[i]);
            //    this->m_dma_hls_fft_run_hdls[i].start();
            //    cols_buff_async_hdls.push_back(buffers_out[i].async("fftColsGraph[" + std::to_string(i) + "].gmio_out", 
            //                                                        XCL_BO_SYNC_BO_AIE_TO_GMIO, 
            //                                                        BLOCK_SIZE_BYTES, 
            //                                                        offset));
            //}


            
            //TODO: THERE SEEMS TO BE A PROBLEM WITH HOW LONG I WAIT (OR DON'T WAIT) HERE FOR THE AIE TO PASS BACK THE DATA UPON EVERY PULSE. 
            //NORMALLY I WOULD JUST USE THE WAIT CALL, BUT MOST OF THE TIME IT SEEMS LIKE IT BLOCKS FOR WAY TOO LONG, AND SOMETIEMS (MORE RARELY)
            //IT DOESN'T BLOCK LONG ENOUGH. I FOUND THAT IT'S BETTER TO EXPERAMENT WITH 0-3 MICROSECONDS. WEIRDLY ENOUGH, SOMETIMES IT FAILS ON
            //0 MICROSECONDS AND IF I BUILD THE APP AGAIN WITHOUT CHANGING ANYTHING, IT WORKS....NOT SURE WHY THIS SEEMS TO BE THE CASE.
            //std::this_thread::sleep_for(std::chrono::microseconds(1));
            //
            //for(int i=0; i<buff_async_hdls.size(); i++) {
            //    buff_async_hdls[i].wait();
            //    //printf("%d: Waiting (buff_size=%d)...\n", i, buff_async_hdls.size());
            //}
            //buff_async_hdls.clear();
        }
        
        for(int sw_id = 0; sw_id < AIE_SWITCHES; sw_id++) {
            int pl_kern_id = buff_idx*m_instances + sw_id;

            // Pass xrt::bo pointer to PL kernel's second arg (ddr_mem) so it can write the AIE data to it
            m_dma_pkt_router_run_hdls[pl_kern_id].set_arg(1, buffers_img_out[buff_idx]);

            // Launch the PL kernel
            m_dma_pkt_router_run_hdls[pl_kern_id].start();
        }

        // Wait for the PL kernels to finish
        //for(int sw_id = 0; sw_id < AIE_SWITCHES; sw_id++) {
        //    m_dma_pkt_router_run_hdls[buff_idx*m_instances + sw_id].wait();
        //}

        // Copy results back to host from PL kernel
        buffers_img_out[buff_idx].sync(XCL_BO_SYNC_BO_FROM_DEVICE, PULSES*RC_SAMPLES*sizeof(cfloat), 0);
    }


    // Block until AIE has finished with above operations. Could do other work on 
    // processor here if needed.
    //for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
    //    //printf("Waiting on idx%d\n", idx);
    //    buff_async_hdls[idx].wait();
    //}
    //std::this_thread::sleep_for(std::chrono::seconds(1));
}

//void SARBackproject::fft(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
//    
//    /* 
//     * The FFT operation uses 2^TP_PARALLEL_POWER subframe processors in parallel each performing 
//     * an FFT of N/(2^TP_PARALLEL_POWER) point size. These subframe outputs are then combined by 
//     * TP_PARALLEL_POWER stages of radix2 to create the final result. As a result, the FFT input
//     * requires the row of data to be split into each port. 
//     *
//     * For example, for 8192 point size and a TP_PARALLEL_POWER = 2, the size for each port
//     * becomes 8192/4=2048. This means: 
//     *
//     * N[0-2047]    -> port0
//     * N[2048-4095] -> port1
//     * N[4096-6143] -> port2 
//     * N[6144-8191] -> port3
//     *
//     * Because of the radix stages, the output is also mixed. Continuing with the example
//     * above, assume the FFT output provides an array of P[8192]. To get the actual order
//     * of the FFT (Q), the array would need to be rearranged as follows:
//     *
//     * P[0]  -> Q[0]
//     * P[4]  -> Q[1]
//     * P[8]  -> Q[2]
//     * P[12] -> Q[3]
//     * .
//     * .
//     * .
//     * P[1]  -> Q[2048]
//     * P[5]  -> Q[2049]
//     * P[9]  -> Q[2050]
//     * P[13] -> Q[2051]
//     * .
//     * .
//     * .
//     * P[2]  -> Q[4096]
//     * P[6]  -> Q[4097]
//     * P[10] -> Q[4098]
//     * P[14] -> Q[4099]
//     * .
//     * .
//     * .
//     * P[3]  -> Q[6144]
//     * P[7]  -> Q[6145]
//     * P[11] -> Q[6146]
//     * P[15] -> Q[6147]
//     *
//     * The reshapeMatrix() function can be used for reshaping the data for the input ports
//     * and strideCols() function can be used for reshaping the data for the output ports.
//     */
//
//    std::vector<xrt::bo::async_handle> buff_async_hdls;
//
//    // 54 MS TO COMPLETE
//    int per_ssr_byte_size = BLOCK_SIZE_BYTES / FFT_NPORTS;
//
//    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
//        for (int ssr = 0; ssr < FFT_NPORTS; ssr++) {
//            buffers_in[buff_idx].async("fftGraph[" + std::to_string(buff_idx) + "].gmio_in[" + std::to_string(ssr) + "]", 
//                                       XCL_BO_SYNC_BO_GMIO_TO_AIE, 
//                                       per_ssr_byte_size, 
//                                       ssr*per_ssr_byte_size);
//
//            buff_async_hdls.push_back(buffers_out[buff_idx].async("fftGraph[" + std::to_string(buff_idx) + "].gmio_out[" + std::to_string(ssr) + "]",
//                                                                  XCL_BO_SYNC_BO_AIE_TO_GMIO, 
//                                                                  per_ssr_byte_size, 
//                                                                  ssr*per_ssr_byte_size));
//        }
//    }
//
//    // Block until AIE has finished with above operations. Could do other work on 
//    // processor here if needed.
//    for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
//        buff_async_hdls[idx].wait();
//    }
//}
//
//void SARBackproject::ifft(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
//
//    std::vector<xrt::bo::async_handle> buff_async_hdls;
//    int per_ssr_byte_size = BLOCK_SIZE_BYTES / FFT_NPORTS;
//
//    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
//        for (int ssr = 0; ssr < FFT_NPORTS; ssr++) {
//            buffers_in[buff_idx].async("ifftGraph[" + std::to_string(buff_idx) + "].gmio_in[" + std::to_string(ssr) + "]", 
//                                       XCL_BO_SYNC_BO_GMIO_TO_AIE, 
//                                       per_ssr_byte_size, 
//                                       ssr*per_ssr_byte_size);
//
//            buff_async_hdls.push_back(buffers_out[buff_idx].async("ifftGraph[" + std::to_string(buff_idx) + "].gmio_out[" + std::to_string(ssr) + "]",
//                                                                  XCL_BO_SYNC_BO_AIE_TO_GMIO, 
//                                                                  per_ssr_byte_size, 
//                                                                  ssr*per_ssr_byte_size));
//        }
//    }
//
//    // Block until AIE has finished with above operations. Could do other work on 
//    // processor here if needed.
//    for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
//        buff_async_hdls[idx].wait();
//    }
//}
//
//void SARBackproject::cplxConj(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers) {
//    int offset = 0;
//    std::vector<xrt::bo::async_handle> buff_async_hdls;
//
//    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
//        buffers_in[buff_idx].async("cplxConjGraph[" + std::to_string(buff_idx) + "].gmio_in", 
//                                   XCL_BO_SYNC_BO_GMIO_TO_AIE, 
//                                   BLOCK_SIZE_BYTES, 
//                                   offset);
//        buff_async_hdls.push_back(buffers_out[buff_idx].async("cplxConjGraph[" + std::to_string(buff_idx) + "].gmio_out",
//                                                              XCL_BO_SYNC_BO_AIE_TO_GMIO, 
//                                                              BLOCK_SIZE_BYTES, 
//                                                              offset));
//    }
//
//    // Block until AIE has finished with above operations. Could do other work on 
//    // processor here if needed.
//    for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
//        buff_async_hdls[idx].wait();
//    }
//}
//
//void SARBackproject::elemMatMult(xrt::aie::bo* buffersA_in, xrt::aie::bo* buffersB_in,
//                                  xrt::aie::bo* buffers_out, int num_of_buffers) {
//    std::vector<xrt::bo::async_handle> buff_async_hdls;
//    int per_ssr_byte_size = BLOCK_SIZE_BYTES / TP_SSR;
//
//    for(int buff_idx = 0; buff_idx < num_of_buffers; buff_idx++) {
//        for (int ssr = 0; ssr < TP_SSR; ssr++) {
//            buffersA_in[buff_idx].async("hpGraph[" + std::to_string(buff_idx) + "].gmio_in_A[" + std::to_string(ssr) + "]", 
//                                        XCL_BO_SYNC_BO_GMIO_TO_AIE, 
//                                        per_ssr_byte_size, 
//                                        ssr*per_ssr_byte_size);
//
//            buffersB_in[buff_idx].async("hpGraph[" + std::to_string(buff_idx) + "].gmio_in_B[" + std::to_string(ssr) + "]", 
//                                        XCL_BO_SYNC_BO_GMIO_TO_AIE, 
//                                        per_ssr_byte_size, 
//                                        ssr*per_ssr_byte_size);
//
//            buff_async_hdls.push_back(buffers_out[buff_idx].async("hpGraph[" + std::to_string(buff_idx) + "].gmio_out[" + std::to_string(ssr) + "]",
//                                                                  XCL_BO_SYNC_BO_AIE_TO_GMIO, 
//                                                                  per_ssr_byte_size, 
//                                                                  ssr*per_ssr_byte_size));
//        }
//    }
//    
//    // Block until AIE has finished with above operations. Could do other work on 
//    // processor here if needed.
//    for(int idx = 0; idx < buff_async_hdls.size(); idx++) {
//        buff_async_hdls[idx].wait();
//    }
//}
