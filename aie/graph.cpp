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
    //std::getline(file, line); // skip header line
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
    //std::getline(file, line); // skip header line
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

    float* xyz_px_array = (float*) GMIO::malloc(PULSES*RC_SAMPLES*sizeof(float)*3);
    int idx = 0;
    for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
        for(int rng_idx = 0; rng_idx < RC_SAMPLES; rng_idx++) {
            
            // X target pixels
            xyz_px_array[idx++] = (rng_idx-HALF_RANGE_SAMPLES)*RANGE_RES;

            // Y target pixels
            xyz_px_array[idx++] = az_res*pulse_idx - half_az_width;

            // Z target pixels
            xyz_px_array[idx++] = 0.0;

            printf("xyz_px_array[%d] = {%f, %f, %f}\n", (idx-3)/3, xyz_px_array[idx-3], xyz_px_array[idx-2], xyz_px_array[idx-1]);
        }
    }

    //// 2) Allocate a second array to hold the correct ordering:
    //float* reordered = (float*)GMIO::malloc(PULSES * RC_SAMPLES * 3 * sizeof(float));
    //
    //// 3) Compute how many blocks per row:
    //const int BLOCK = 16;
    //const int NUM_BLOCKS = RC_SAMPLES / BLOCK;  // e.g. 64/16 = 4
    //
    //// 4) For every row j and every block b, memcpy that 16-point chunk into its new slot:
    //for (int j = 0; j < PULSES; j++) {
    //    for (int b = 0; b < NUM_BLOCKS; b++) {
    //        // Original chunk starts at row j, block b of size BLOCK points:
    //        size_t original_offset_bytes = ((size_t)j * RC_SAMPLES + (size_t)b * BLOCK) * 3 * sizeof(float);
    //    
    //        // In the reordered array, this should land at:
    //        size_t target_offset_bytes = ((size_t)b * PULSES * BLOCK + (size_t)j * BLOCK) * 3 * sizeof(float);
    //    
    //        // Copy BLOCK points (each point = 3 floats):
    //        memcpy(reinterpret_cast<char*>(reordered) + target_offset_bytes,
    //               reinterpret_cast<char*>(xyz_px_array) + original_offset_bytes,
    //               (size_t)BLOCK * 3 * sizeof(float));
    //    }
    //}
    //
    //// 5) (Optional) Verify a few entries:
    //for (int i = 0; i < PULSES * RC_SAMPLES; i++) {
    //    float x = reordered[i * 3 + 0];
    //    float y = reordered[i * 3 + 1];
    //    float z = reordered[i * 3 + 2];
    //    printf("reordered[%d] = {%f, %f, %f}\n", i, x, y, z);
    //}

    //int px_per_chunk = 16;
    //int total_px = PULSES*RC_SAMPLES;
    //int px_per_solver = total_px/IMG_SOLVERS;
    //int chunks_per_solver = px_per_solver/px_per_chunk;
    //int chunks_per_switch = IMG_SOLVERS_PER_SWITCH*chunks_per_solver;
    //int chunks_per_pulse = RC_SAMPLES/px_per_chunk;

    //float* xyz_px_array = (float*) GMIO::malloc(total_px*sizeof(float)*3);

    //// Create a 3D view of xyz_px_array so we can index [solver][pixel_in_solver][coord]
    ////float (*xyz_px_array_3d)[px_per_solver][3] = (float (*)[px_per_solver][3]) xyz_px_array;

    //int global_px = 0;
    //for(int range_block = 0; range_block < chunks_per_pulse; range_block++) {
    //    for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
    //        for(int px_in_chunk = 0; px_in_chunk < px_per_chunk; px_in_chunk++) {
    //        
    //            // 1) Derive X, Y, and Z target pixels
    //            int rng_idx = range_block * px_per_chunk + px_in_chunk;
    //            float x_px = (rng_idx-HALF_RANGE_SAMPLES)*RANGE_RES;
    //            float y_px = az_res*pulse_idx - half_az_width;
    //            float z_px = 0.0;


    //            // 2) Which chunk the global_px belongs to
    //            int chunk_num = global_px / px_per_chunk;

    //            // 3) Which switch cluster the chunk is in
    //            int switch_idx = chunk_num / chunks_per_switch;

    //            // 4) Within this switch, which chunk-slot is it
    //            int within_switch = chunk_num % chunks_per_switch;

    //            // 5) Which solver within that switch
    //            int solver_within_switch = within_switch % IMG_SOLVERS_PER_SWITCH;

    //            // 6) Which chunk-index inside that solver
    //            int chunk_within_solver = within_switch / IMG_SOLVERS_PER_SWITCH;  

    //            // 7) Absolute solver id
    //            int solver_id = switch_idx * IMG_SOLVERS_PER_SWITCH + solver_within_switch;
    //    
    //            // 8) Within that solver's block, which pixel index
    //            //int px_in_chunk = global_px % px_per_chunk;
    //            int px_idx_in_solver = chunk_within_solver * px_per_chunk + px_in_chunk;



    //            // 9) Sort in X, Y and Z 3D array
    //            int xyz_idx = solver_id * px_per_solver + px_idx_in_solver;
    //            int float_xyz_idx = xyz_idx*3;
    //            //xyz_px_array_3d[solver_id][px_idx_in_solver][0] = x_px;
    //            //xyz_px_array_3d[solver_id][px_idx_in_solver][1] = y_px;
    //            //xyz_px_array_3d[solver_id][px_idx_in_solver][2] = z_px;
    //            xyz_px_array[float_xyz_idx+0] = x_px;
    //            xyz_px_array[float_xyz_idx+1] = y_px;
    //            xyz_px_array[float_xyz_idx+2] = z_px;

    //            // 10) Advance global_pixel
    //            global_px++;

    //            printf("range_blk: %d | pulse_idx: %d | px_in_chunk: %d | chunk_num: %d | switch_idx: %d | within_switch: %d | solver_within_switch: %d | chunk_within_solver: %d | solver_id: %d | px_idx_in_solver: %d | xyz_idx: %d | x: %f | y: %f\n", 
    //                    range_block,
    //                    pulse_idx,
    //                    px_in_chunk,
    //                    chunk_num, 
    //                    switch_idx, 
    //                    within_switch, 
    //                    solver_within_switch, 
    //                    chunk_within_solver,
    //                    solver_id,
    //                    px_idx_in_solver,
    //                    xyz_idx,
    //                    x_px,
    //                    y_px);
    //        }

    //    }
    //}

    //for (int px = 0;  px < total_px;  px++) {
    //    int idx = px * 3;
    //    float x = xyz_px_array[idx + 0];
    //    float y = xyz_px_array[idx + 1];
    //    float z = xyz_px_array[idx + 2];
    //    printf("xyz_px_array[%d] = {%f, %f, %f}\n", idx/3, x, y, z);
    //}
    
    //int idx = 0;
    //for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
    //    for(int rng_idx = 0; rng_idx < RC_SAMPLES; rng_idx++) {
    //        idx+=3;
    //        printf("xyz_px_array[%d] = {%f, %f, %f}\n", (idx-3)/3, xyz_px_array[idx-3], xyz_px_array[idx-2], xyz_px_array[idx-1]);
    //    }
    //}


    // Number of iteration for the AIE graphs to run
    for (int inst = 0; inst < INSTANCES; inst++) {
        bpGraph[inst].run(PULSES);
    }

    // Number of target pixels for each px_demux_kern to process
    //int px_per_demux_kern = ((PULSES*RC_SAMPLES)/AIE_SWITCHES);

    // Loop through pipeline ITER times
    int inst = 0;
    
    for(int iter=0; iter<ITER; iter++) {
        printf("\nPERFORM BACKPROJECTION (ITER = %d) (INST = %d)\n", iter, inst);

        // Pass in slowtime data into AI kernels
        bpGraph[inst].gmio_in_st.gm2aie_nb(broadcast_data_array, PULSES*BC_ELEMENTS*sizeof(float));

        // Pass in other data into bp AI kernels
        for(int pulse_idx=0; pulse_idx<PULSES; pulse_idx++) {

            bpGraph[inst].gmio_in_rc.gm2aie_nb(rc_array + pulse_idx*RC_SAMPLES, RC_SAMPLES*sizeof(cfloat));
            //for (int sw_id=0; sw_id<AIE_SWITCHES; sw_id++) {
            //    bpGraph[inst].bpCluster[sw_id].gmio_in_xyz_px.gm2aie_nb(xyz_px_array + sw_id*px_per_demux_kern*3, px_per_demux_kern*sizeof(float)*3);
            //}

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
