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
, m_st_dataset_filename(st_dataset_filename)
, m_rc_dataset_filename(rc_dataset_filename)
, m_img_out_filename(img_out_filename)
, m_iter(iter)
, m_instances(instances)
, m_device(0)
, m_uuid(m_device.load_xclbin(this->m_xclbin_filename))
, m_broadcast_data_buffer(m_device, PULSES*BC_ELEMENTS*sizeof(float), xrt::bo::flags::normal, 0)
, m_broadcast_data_array(m_broadcast_data_buffer.map<float*>())
, m_rc_buffer(m_device, PULSES*RC_SAMPLES*sizeof(cfloat), xrt::bo::flags::normal, 0)
, m_rc_array(m_rc_buffer.map<cfloat*>())
{
    // Instantiate kernels, handlers, and buffers for multiple instances
    // TODO: Instances wont work right if it is greater than 1. This is fine
    // for now since there are not plans to instantiate the design more than
    // once, but would be nice if this still worked for future proofing
    for(int i = 0; i < this->m_instances; i++) {
        std::string dma_stride_kernel_str = "dma_stride_controller:{dma_stride_controller_" + std::to_string(i) + "}";
        m_dma_stride_kernels.push_back(xrt::kernel(m_device, m_uuid, dma_stride_kernel_str));
        m_dma_stride_run_hdls.push_back(xrt::run(m_dma_stride_kernels[i]));

        for(int sw_id = 0; sw_id < AIE_SWITCHES; sw_id++) {
            std::string dma_pkt_router_kernel_str = "dma_pkt_router:{dma_pkt_router_" + std::to_string(sw_id) + "}";
            m_dma_pkt_router_kernels.push_back(xrt::kernel(m_device, m_uuid, dma_pkt_router_kernel_str));
            m_dma_pkt_router_run_hdls.push_back(xrt::run(m_dma_pkt_router_kernels[sw_id]));
        }
        
        m_xyz_px_buffers.push_back(xrt::bo(m_device, PULSES*AZ_SAMPLES*RC_SAMPLES*2*sizeof(uint64_t), m_dma_stride_kernels[0].group_id(0)));
        m_xyz_px_arrays.push_back(m_xyz_px_buffers[i].map<uint64_t*>());
        
        // All kernel instances map to the same DDR memory bank. This means
        // kernel1.group_id(1) should return the same bank ID as
        // kernel0.group_id(1). Therefore, it doesnt matter which instances
        // group_id is used, the result is the same memory bank.
        m_img_buffers.push_back(xrt::bo(m_device, AZ_SAMPLES*RC_SAMPLES*8, m_dma_pkt_router_kernels[0].group_id(1)));
        m_img_arrays.push_back(m_img_buffers[i].map<cfloat*>());

        std::string bpGraph_str = "bpGraph[" + std::to_string(i) + "]";
        m_bp_graph_hdls.push_back(xrt::graph(this->m_device, this->m_uuid, bpGraph_str));
    }
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
    double delta_az;
    double total_az;
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
        delta_az = std::fabs(mean_diff);
        double min_az = *std::min_element(az_ant, az_ant + PULSES);
        double max_az = *std::max_element(az_ant, az_ant + PULSES);
        total_az = max_az - min_az;
        az_res = C/(2.0*total_az*MIN_FREQ);
        double az_width = C/(2.0*delta_az*MIN_FREQ);
        half_az_width = az_width/2.0;
    }

    // Determine the maximum scene size of the image (m)
    double max_wr = C/(2*RANGE_FREQ_STEP);
    double max_wx = C/(2*delta_az*MIN_FREQ);

    // Determine the resolution of the image (m)
    double dr = C/(2*RANGE_FREQ_STEP*424);
    double dx = C/(2*total_az*MIN_FREQ);

    printf("Total Az: %.1f deg, %.4f rad\n", total_az*(180.0/PI), total_az);
    printf("Maximum Scene Size: %.2f m range, %.2f m cross-range\n", max_wr, max_wx);
    printf("Resolution: %.2fm range, %.2f m cross-range\n", dr, dx);

    // Pack each pair of 32-bit floats into one 64-bit word
    int word_idx = 0;
    for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
        for(int az_idx = 0; az_idx < AZ_SAMPLES; az_idx++) {
            for(int rng_idx = 0; rng_idx < RC_SAMPLES; rng_idx++) {
                
                // X target pixels
                float x = (rng_idx-HALF_RANGE_SAMPLES)*RANGE_RES;

                // Y target pixels
                float y = az_res*az_idx - half_az_width;

                // Z target pixels
                float z = 0.0;

                // **** PACK X AND Y INTO WORD0 **** //
                uint32_t bits_x, bits_y, bits_z;
                memcpy(&bits_x, &x, sizeof(float));
                memcpy(&bits_y, &y, sizeof(float));

                // Lower 32 bits = x | Higher 32 bits = y
                uint64_t w0 = ((uint64_t)bits_y << 32) | bits_x;
                m_xyz_px_arrays[0][word_idx++] = w0;

                // **** PACK Z AND 0 INTO WORD1 **** //
                memcpy(&bits_z, &z, sizeof(float));

                // Lower 32 bits = z | Higher 32 bits = 0
                uint64_t w1 = bits_z;
                m_xyz_px_arrays[0][word_idx++] = w1;
            }
        }
    }

    m_xyz_px_buffers[0].sync(XCL_BO_SYNC_BO_TO_DEVICE, PULSES*AZ_SAMPLES*RC_SAMPLES*2*sizeof(uint64_t), 0);

    // Pass xrt::bo pointer to PL kernel's first arg (ddr_mem) so it can read the pixel data from it
    m_dma_stride_run_hdls[0].set_arg(0, m_xyz_px_buffers[0]);

    // Launch the DMA Stride Controller PL kernel
    m_dma_stride_run_hdls[0].start();
}

void SARBackproject::runGraphs() {
    // Run graphs indefinitley (instead of for specific amount of iterations)
    for(int i = 0; i < this->m_instances; i++) {
        this->m_bp_graph_hdls[i].run(0);
    }
}

void SARBackproject::bp(xrt::aie::bo* buffers_broadcast_data_in, xrt::aie::bo* buffers_rc_in, 
        xrt::bo* buffers_img_out, int num_of_buffers) {
    
    std::vector<xrt::bo::async_handle> buff_async_hdls;
    
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
            

            for (int kern_id = 0; kern_id < IMG_SOLVERS; kern_id++) {
                // Dump image if on last pulse, otherwise keep focusing the image
                if (pulse_idx == PULSES-1) {
                    this->m_bp_graph_hdls[buff_idx].update("bpGraph[" + std::to_string(buff_idx) + "].rtp_dump_img_in[" + std::to_string(kern_id) + "]", 1);
                } else {
                    this->m_bp_graph_hdls[buff_idx].update("bpGraph[" + std::to_string(buff_idx) + "].rtp_dump_img_in[" + std::to_string(kern_id) + "]", 0);
                }
            }
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
        buffers_img_out[buff_idx].sync(XCL_BO_SYNC_BO_FROM_DEVICE, AZ_SAMPLES*RC_SAMPLES*sizeof(cfloat), 0);
    }
}

