// By: Austin Owens
// Date: 5/16/2025
// Desc: Performs SAR backprojection

#include "graph.h"

uint8_t bp_graph_insts = 0;

const int INSTANCES = 1;
BackProjectionGraph bpGraph[INSTANCES];

#if defined(__AIESIM__) || defined(__X86SIM__)
#include <chrono>
#include <thread>
#include <regex>
#include <unistd.h>

void unwrap(double* angles) {
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

int main(int argc, char ** argv) {

    // Number of iterations to run
    const int ITER = 1;

    // Generate random seed
    std::time_t seed = std::time(0);
    std::srand(seed);
    std::cout << "Rand Seed: " << seed << std::endl;

    // Total error count
    int total_err_cnt = 0;

    // Initialize all AIE graphs
    for(int inst=0; inst<INSTANCES; inst++) {
        bpGraph[inst].init();
    }

    // OPEN SAR DATASET FILES
    float* broadcast_data_array = (float*) GMIO::malloc(PULSES*BC_ELEMENTS*sizeof(float));
    
    // Current working dir inside build/hw/aiesim/
    std::ifstream st_file("../../../design/test_data/gotcha_slowtime_pass1_360deg_HH.csv");
    if (!st_file.is_open()) {
        std::cerr << "Error opening st_file dataset!" << std::endl;
        return 1;
    }
    std::string line;
    int pulse_idx = 0;
    while (std::getline(st_file, line) && pulse_idx < PULSES) {
        std::stringstream ss(line);
        std::string value;

        // Read each column from the line and store in respective array
        std::getline(ss, value, ',');
        broadcast_data_array[BC_ELEMENTS*pulse_idx] = std::stof(value);

        std::getline(ss, value, ',');
        broadcast_data_array[BC_ELEMENTS*pulse_idx + 1] = std::stof(value);

        std::getline(ss, value, ',');
        broadcast_data_array[BC_ELEMENTS*pulse_idx + 2] = std::stof(value);

        std::getline(ss, value, ',');
        broadcast_data_array[BC_ELEMENTS*pulse_idx + 3] = std::stof(value);

        pulse_idx++;
    }

    cfloat* rc_array = (cfloat*) GMIO::malloc(PULSES*RC_SAMPLES*sizeof(cfloat));
     
    // Current working dir inside build/hw/aiesim/
    std::string rc_filename = "../../../design/test_data/gotcha_" + 
                              std::to_string(RC_SAMPLES) + 
                              "-out-of-424-rc-samples_pass1_360deg_HH.csv";
    std::ifstream rc_file(rc_filename);
    if (!rc_file.is_open()) {
        std::cerr << "Error opening rc_file dataset!" << std::endl;
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
            rc_array[pulse_idx*RC_SAMPLES + rc_samp_cnt] = (cfloat) {real_part, imag_part};

            rc_samp_cnt++;
        }
        rc_samp_cnt = 0;
        pulse_idx++;
    }

    // AZIMUTH RESOLUTION GRID
    double az_res = 0;
    double half_az_width = 0;
    if (PULSES != 1) {
        double az_ant[PULSES];
        for (int i = 0; i < PULSES; i++) {
            // atan2(y_ant_pos / x_ant_pos)
            az_ant[i] = std::atan2(broadcast_data_array[BC_ELEMENTS*i+1], broadcast_data_array[BC_ELEMENTS*i]);
        }
        unwrap(az_ant);
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

    int idx = 0;
    float* tmp_px_array = (float*) malloc(PULSES*RC_SAMPLES*sizeof(float)*3);
    for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
        for(int rng_idx = 0; rng_idx < RC_SAMPLES; rng_idx++) {
            
            // X target pixels
            tmp_px_array[idx++] = (rng_idx-HALF_RANGE_SAMPLES)*RANGE_RES;

            // Y target pixels
            tmp_px_array[idx++] = az_res*pulse_idx - half_az_width;

            // Z target pixels
            tmp_px_array[idx++] = 0.0;

            printf("tmp_px_array[%d] = {%f, %f, %f}\n", (idx-3)/3, tmp_px_array[idx-3], tmp_px_array[idx-2], tmp_px_array[idx-1]);
        }
    }

    // Number of target pixels in a single pulse
    const int px_per_pulse = AZ_SAMPLES * RC_SAMPLES;

    // Number of target pixels destined for each switch per pulse
    const int px_per_switch = px_per_pulse/AIE_SWITCHES;

    // Number of target pixels destined for each AIE kernel per pulse
    const int px_per_kern = px_per_pulse/IMG_SOLVERS;
    
    float* xyz_px_array = (float*) GMIO::malloc(PULSES*RC_SAMPLES*sizeof(float)*3);
    int idx_px_component;
    int new_idx = 0;

    // Iterate over every pulse
    //for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
    for(int pulse_idx = 0; pulse_idx < 1; pulse_idx++) {

        // Iterate over each switch
        for(int sw_num = 0; sw_num < AIE_SWITCHES; sw_num++) {

            // Starting pixel index on each pulse
            int start_px_idx_per_pulse = sw_num*px_per_switch + pulse_idx*px_per_pulse;

            // Stride 16 pixels
            for(int base_px = start_px_idx_per_pulse; base_px < px_per_kern + start_px_idx_per_pulse; base_px += 16) {

                // Strides the number of pixels that will be sent to an AIE kernel
                for(int start_px = base_px; start_px < px_per_switch + start_px_idx_per_pulse; start_px += px_per_kern) {

                    // Extract 16 pixels
                    for(int idx_px = start_px; idx_px < start_px + 16; idx_px++) {

                        // Because idx_128b is a 128 bit sample index representing
                        // a pixel, but ddr_mem is 64bit, we need to multiply by 2
                        idx_px_component = idx_px*3;

                        xyz_px_array[new_idx++] = tmp_px_array[idx_px_component];
                        xyz_px_array[new_idx++] = tmp_px_array[idx_px_component+1];
                        xyz_px_array[new_idx++] = tmp_px_array[idx_px_component+2];

                        printf("orig_idx: %d | xyz_px_array[%d] = {%f, %f, %f}\n", 
                                idx_px, new_idx-3, xyz_px_array[new_idx-3], xyz_px_array[new_idx-2], xyz_px_array[new_idx-1]);
                        
                    }
                }
            }
        }
    }
    free(tmp_px_array);


    // Number of iteration for the AIE graphs to run
    for (int inst = 0; inst < INSTANCES; inst++) {
        bpGraph[inst].run(PULSES);
    }

    // Number of target pixels for each px_demux_kern to process
    int px_per_demux_kern = ((PULSES*RC_SAMPLES)/AIE_SWITCHES);

    // Loop through pipeline ITER times
    int inst = 0;
    
    for(int iter=0; iter<ITER; iter++) {
        printf("\nPERFORM BACKPROJECTION (ITER = %d) (INST = %d)\n", iter, inst);

        // Pass in slowtime data into AI kernels
        bpGraph[inst].gmio_in_st.gm2aie_nb(broadcast_data_array, PULSES*BC_ELEMENTS*sizeof(float));

        // Pass in other data into bp AI kernels
        for(int pulse_idx=0; pulse_idx<PULSES; pulse_idx++) {

            bpGraph[inst].gmio_in_rc.gm2aie_nb(rc_array + pulse_idx*RC_SAMPLES, RC_SAMPLES*sizeof(cfloat));
            for (int sw_id=0; sw_id<AIE_SWITCHES; sw_id++) {
                bpGraph[inst].bpCluster[sw_id].gmio_in_xyz_px.gm2aie_nb(xyz_px_array + sw_id*px_per_demux_kern*3, px_per_demux_kern*sizeof(float)*3);
            }

            for(int kern_id=0; kern_id<IMG_SOLVERS; kern_id++) {
                // Dump image if on last pulse, otherwise keep focusing the image
                if (pulse_idx == PULSES-1)
                    bpGraph[inst].update(bpGraph[inst].rtp_dump_img_in[kern_id], 1);
                else
                    bpGraph[inst].update(bpGraph[inst].rtp_dump_img_in[kern_id], 0);
            }
        }
    }

    // Block is needed for PLIO output to completley write to file (otherwise
    // it will be blank) 
    bpGraph[inst].wait();
    bpGraph[inst].end();
    return total_err_cnt;
}

#endif
