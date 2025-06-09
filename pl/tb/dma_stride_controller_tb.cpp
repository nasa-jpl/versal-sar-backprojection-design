// By: Austin Owens
// Date: 5/24/2025
// Desc: Test DMA stride controller Programmable Logic (PL) kernel. The DMA
// stride controller passes data into the AIE in a specific order to better
// accomidate the AIE kernels operating on the data.

#include "../dma_stride_controller.h"
#include <fstream>
#include <iostream>
#include <algorithm>

using namespace std;

// It helps to imagine this testbench from the AI Engine perspective

void unwrap(double* angles) {
    // Store the first original angle for difference calculations.
    double prev_orig = angles[0];

    // The first angle remains unchanged.
    for (int i = 1; i < AZ_SAMPLES; i++) {

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

void gen_grid_params(double &az_res, double &half_az_width) {
    // OPEN SAR DATASET FILES
    float* broadcast_data_array = (float*) malloc(PULSES*BC_ELEMENTS*sizeof(float));
    
    // Current working dir inside versal-design-build/build/hw/plsim/dma_pkt_router_testbench/solution1/csim/build
    ifstream st_file("../../../../../../../design/test_data/gotcha_slowtime_pass1_360deg_HH.csv");
    if (!st_file.is_open()) {
        cerr << "Error opening st_file dataset!" << endl;
        exit(-1);
    }
    string line;
    int pulse_idx = 0;
    //getline(file, line); // skip header line
    while (getline(st_file, line) && pulse_idx < PULSES) {
        stringstream ss(line);
        string value;

        // Read each column from the line and store in respective array
        getline(ss, value, ',');
        broadcast_data_array[BC_ELEMENTS*pulse_idx] = stof(value);

        getline(ss, value, ',');
        broadcast_data_array[BC_ELEMENTS*pulse_idx + 1] = stof(value);

        getline(ss, value, ',');
        broadcast_data_array[BC_ELEMENTS*pulse_idx + 2] = stof(value);

        getline(ss, value, ',');
        broadcast_data_array[BC_ELEMENTS*pulse_idx + 3] = stof(value);

        pulse_idx++;
    }

    // AZIMUTH RESOLUTION GRID
    if (PULSES != 1) {
        double az_ant[PULSES];
        for (int i = 0; i < PULSES; i++) {
            // atan2(y_ant_pos / x_ant_pos)
            az_ant[i] = atan2(broadcast_data_array[BC_ELEMENTS*i+1], broadcast_data_array[BC_ELEMENTS*i]);
        }
        unwrap(az_ant);
        double sum_diff = 0.0;
        for (int i = 1; i < PULSES; i++) {
            sum_diff += (az_ant[i] - az_ant[i - 1]);
        }
        double mean_diff = sum_diff / (PULSES - 1);
        double delta_az = fabs(mean_diff);
        double min_az = *min_element(az_ant, az_ant + PULSES);
        double max_az = *max_element(az_ant, az_ant + PULSES);
        double total_az = max_az - min_az;
        az_res = C/(2.0*total_az*MIN_FREQ);
        double az_width = C/(2.0*delta_az*MIN_FREQ);
        half_az_width = az_width/2.0;
    }
    free(broadcast_data_array);

}

int main() {
    // Total number of AXIS 128b samples
    const int AXIS128_SAMPLES = (PULSES*AZ_SAMPLES*RC_SAMPLES)/AIE_SWITCHES;
    
    // DDR is a 64bit memory space. Because the pixels we want to store consist
    // of X, Y, and Z floats, that takes up 12 bytes per sample (4 bytes per
    // component). We can take better advantage of the space in the DDR by
    // having X and Y components take up the lower 32 bits and upper 32 bits of
    // a 64 bit memory space, then have the Z component take up the next 64 bit
    // memory space (while padding the higher 32 bits). Therefore, we only need
    // PULSES*AZ_SAMPLES*RC_SAMPLES*2 UINT64 memory spaces instead of
    // PULSES*AZ_SAMPLES*RC_SAMPLES*3.
    const int DDR_UINT64_SAMPLES = PULSES*AZ_SAMPLES*RC_SAMPLES*2;

    // Allocate DDR memory buffer
    ap_uint<64>* ddr_mem = (ap_uint<64>*) malloc(DDR_UINT64_SAMPLES * sizeof(ap_uint<64>));
    
    // Get grid params
    double az_res = 0;
    double half_az_width = 0;
    gen_grid_params(az_res, half_az_width);

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
                ap_uint<64> w0 = 0;
                w0.range(31, 0) = bits_x;
                w0.range(63, 32) = bits_y;
                ddr_mem[word_idx++] = w0;

                // **** PACK Z AND 0 INTO WORD1 **** //
                memcpy(&bits_z, &z, sizeof(float));

                // Lower 32 bits = z | Higher 32 bits = 0
                ap_uint<64> w1 = 0;
                w1.range(31,  0) = bits_z;
                w1.range(63,  32) = 0;
                ddr_mem[word_idx++] = w1;
            }
        }
    }
    

    // Debug print what is in ddr_mem
    for (int i = 0; i < DDR_UINT64_SAMPLES; i++) {
        uint32_t lo =  (uint32_t) ddr_mem[i].range(31, 0);
        uint32_t hi =  (uint32_t) ddr_mem[i].range(63, 32);

        float f_lo, f_hi;
        memcpy(&f_lo, &lo, sizeof(float));
        memcpy(&f_hi, &hi, sizeof(float));

        printf("ddr_mem[%zu] = {%.6f, %.6f}\n", i, f_lo, f_hi);
    }
    
    // Execute PL kernel
    hls::stream<ap_axiu<128, 0, 0, 0>> pl_stream_out[AIE_SWITCHES];
    int ret = dma_stride_controller(ddr_mem, pl_stream_out);

    for(int switch_num=0; switch_num<AIE_SWITCHES; switch_num++) {

        // Current working dir inside versal-design-build/build/hw/plsim/dma_stride_controller_testbench/solution1/csim/build
        // Assumes there is only one instance of bpGraph (hence the 0 on plio_to_aie_switch_0_*)
        string aie_data_str = "../../../../plsimulator_output/plio_to_aie_switch_0_" + to_string(switch_num) + ".csv";
        FILE *trgt_px_fp = fopen(aie_data_str.c_str(), "w");
        if (trgt_px_fp == NULL) {
            perror("Error opening output CSV file");
            return 1;
        }

        fprintf(trgt_px_fp, "CMD, D, D, D, D, TLAST, TKEEP\n");

        //int ddr_offset = DDR_UINT64_SAMPLES_PER_PULSE*(pulse_idx + switch_num*PULSES);
        //printf("ddr_offset: %d\n", ddr_offset);


        ap_axiu<128, 0, 0, 0> pkt;
        int cnt = 0;
        while (!pl_stream_out[switch_num].empty()) {
            // Read one 128-bit AXI packet from the stream
            pkt = pl_stream_out[switch_num].read();

            uint32_t seg0 = (uint32_t) pkt.data.range(31,  0);
            uint32_t seg1 = (uint32_t) pkt.data.range(63, 32);
            uint32_t seg2 = (uint32_t) pkt.data.range(95, 64);
            uint32_t seg3 = (uint32_t) pkt.data.range(127,96);
            int tlast = (int)pkt.last;
            int tkeep = (int)pkt.keep;

            // Reinterpret each 32-bit segment as a float
            float f0, f1, f2, f3;
            memcpy(&f0, &seg0, sizeof(float));
            memcpy(&f1, &seg1, sizeof(float));
            memcpy(&f2, &seg2, sizeof(float));
            memcpy(&f3, &seg3, sizeof(float));

            printf("pl_stream[%d][%d/%d] = {%.6f, %.6f, %.6f, %.6f}, TLAST=%d\n", 
                    switch_num, cnt, AXIS128_SAMPLES-1, f0, f1, f2, f3, tlast);

            fprintf(trgt_px_fp, "DATA, %f, %f, %f, %f, %d, %d\n", f0, f1, f2, f3, tlast, tkeep);

            cnt++;
        }
        fclose(trgt_px_fp);
    }

    free(ddr_mem);

    return 0;
}
