// By: Austin Owens
// Date: 6/3/2024
// Desc: Performs an image cross correlation using FFT's
//
// Notes: 
// This cross correlation implementation solves for a 2D FFT of the image by
// computing the 1D FFT of the rows in the AIE first, followed by a transpose 
// on the host processor, then another 1D FFT on the cols in the AIE. Typically
// when you get the col-wise result back, you would perform another transpose
// to get the completed 2D FFT; I DON'T DO THIS TRANSPOSE. Because I eventually
// need to do a 2D IFFT to get the finalized cross correlated image, I can save
// on two transpose operations by not performing the transpose after receiving
// the col-wise 1D FFT and IFFT responses from the AIE. 

#include "graph.h"

uint8_t fft_graph_insts = 0;
uint8_t ifft_graph_insts = 0;
uint8_t cplx_conj_graph_insts = 0;
uint8_t hp_graph_insts = 0;
uint8_t bp_graph_insts = 0;

const int INSTANCES = 1;
////FFTGraph<2,5> fftGraph[INSTANCES];
////IFFTGraph<22,5> ifftGraph[INSTANCES];
//FFTGraph<-1,-1> fftGraph[INSTANCES];
//IFFTGraph<-1,-1> ifftGraph[INSTANCES];
////CplxConjGraph cplxConjGraph[INSTANCES];
//HPGraph hpGraph[INSTANCES];
BackProjectionGraph bpGraph[INSTANCES];

#if defined(__AIESIM__) || defined(__X86SIM__)
//#include "math.h"
#include <chrono>
#include <thread>
#include <regex>
#include <unistd.h>

//const float EPSILON = 0.0001;
//const int BLOCK_SIZE_ENTRIES = MAT_ROWS*MAT_COLS;
//const int BLOCK_SIZE_BYTES = BLOCK_SIZE_ENTRIES * sizeof(cfloat);

// FFT output should all be FFT_SAMPLE_DATA for both real and imaginary
//int fftErrorCheck(cfloat* data_array) {
//    int error_cnt = 0;
//    float real, imag;
//
//    for(int r = 0; r < MAT_ROWS; r++) {
//        for(int c = 0; c < MAT_COLS; c++) {
//            int index = (r*MAT_COLS)+c;
//            real = data_array[index].real;
//            imag = data_array[index].imag;
//            //printf("FFT[%d, %d] = {%f, %f}\n", r, c, real, imag);
//            if(real != FFT_SAMPLE_DATA || imag != FFT_SAMPLE_DATA) {
//                printf("ERROR FFT[%d, %d] = {%f, %f}\n", r, c, real, imag);
//                error_cnt++;
//            }
//        }
//    }
//
//    printf("FFT ERRORS: %d\n", error_cnt);
//
//    return error_cnt;
//}

//int ifftErrorCheck(cfloat* data_array, cfloat* gold_data_array) {
//    int error_cnt = 0;
//    float real, imag;
//    float gold_real, gold_imag;
//
//    for(int r = 0; r < MAT_ROWS; r++) {
//        for(int c = 0; c < MAT_COLS; c++) {
//            int index = (r*MAT_COLS)+c;
//            real = data_array[index].real/TP_POINT_SIZE;
//            imag = data_array[index].imag/TP_POINT_SIZE;
//            gold_real = gold_data_array[index].real;
//            //gold_imag = -gold_data_array[index].imag; // Adding negative to accommodate for complx conjugation
//            gold_imag = gold_data_array[index].imag; // Adding negative to accommodate for complx conjugation
//            //printf("IFFT[%d, %d] = {%f, %f}\n", r, c, real, imag);
//
//            if(std::abs(real-gold_real) > EPSILON || std::abs(imag-gold_imag) > EPSILON) {
//                printf("ERROR IFFT[%d, %d] = {%f, %f} | EXPECTED {%f, %f}\n", r, c, real, imag, gold_real, gold_imag);
//                error_cnt++;
//            }
//        }
//    }
//
//    printf("IFFT ERRORS: %d\n\n", error_cnt);
//
//    return error_cnt;
//}
//
//void print_arr(cfloat* arr, int num_rows, int num_cols) {
//    for(int r=0; r<num_rows; r++) {
//        for(int c=0; c<num_cols; c++) {
//            int index = (r*MAT_COLS) + c;
//            printf("array[%d][%d] = {%f, %f}\n", r, c, arr[index].real, arr[index].imag);
//        }
//    }
//}
//
//void reorderDataArray(cfloat* data_array) {
//    // Create a temporary array to hold the rearranged data
//    cfloat* temp = (cfloat*) GMIO::malloc(BLOCK_SIZE_BYTES);
//    
//    // Rearrange the data into the temp array
//    int chunk_size = MAT_COLS;
//    int num_datasets = BLOCK_SIZE_ENTRIES / chunk_size;
//
//    for (int dataset = 0; dataset < num_datasets; dataset++) {
//        int offset = dataset * chunk_size;
//
//        for (int i = 0; i < chunk_size; i++) {
//            int new_index = offset + (i / 4) + (i % 4) * (chunk_size / 4);
//            temp[new_index] = data_array[offset + i];
//        }
//    }
//
//    // Copy the rearranged data back into the original array
//    for (int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
//        data_array[i] = temp[i];
//    }
//
//    printf("\nFLIPPED DATA_ARRAY\n");
//    for(int r=0; r<MAT_ROWS; r++) {
//        for(int c=0; c<MAT_COLS; c++) {
//            int index = (r*MAT_COLS) + c;
//            printf("data_array[%d][%d] = {%f, %f}\n", r, c, data_array[index].real, data_array[index].imag);
//        }
//    }
//
//    // Free the temporary array
//    free(temp);
//}

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
        //fftGraph[inst].init();
        //hpGraph[inst].init();
        //ifftGraph[inst].init();
        ////cplxConjGraph[inst].init();
        bpGraph[inst].init();
    }
   
    // Allocating memory
    //std::vector<cfloat*> data_arrays;
    //for(int inst=0; inst<INSTANCES; inst++) {
    //    data_arrays.push_back((cfloat*) GMIO::malloc(BLOCK_SIZE_BYTES));
    //}
    //cfloat* one_array = (cfloat*) GMIO::malloc(BLOCK_SIZE_BYTES);

    // Allocate and populate memory
    //cfloat* gold_data_array = (cfloat*) GMIO::malloc(BLOCK_SIZE_BYTES);
    //for(int r = 0; r < MAT_ROWS; r++) {
    //    for(int c = 0; c < MAT_COLS; c++) {
    //        int index = (r*MAT_COLS)+c;
    //        if (c == 1) {
    //            gold_data_array[index] = (cfloat) {FFT_SAMPLE_DATA, FFT_SAMPLE_DATA};
    //        }
    //        else {
    //            gold_data_array[index] = (cfloat) {0, 0};
    //        }
    //    }
    //}

    // OPEN SAR DATASET FILES
    float* x_ant_data_array = (float*) GMIO::malloc(PULSES*sizeof(float));
    float* y_ant_data_array = (float*) GMIO::malloc(PULSES*sizeof(float));
    float* z_ant_data_array = (float*) GMIO::malloc(PULSES*sizeof(float));
    float* ref_range_data_array = (float*) GMIO::malloc(PULSES*sizeof(float));
    
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
        x_ant_data_array[pulse_idx] = std::stof(value);

        std::getline(ss, value, ',');
        y_ant_data_array[pulse_idx] = std::stof(value);

        std::getline(ss, value, ',');
        z_ant_data_array[pulse_idx] = std::stof(value);

        std::getline(ss, value, ',');
        ref_range_data_array[pulse_idx] = std::stof(value);

        pulse_idx++;
    }
    
    //for(int i=0; i<2; i++) {
    //    printf("x_ant_data_array[%d] = %f\n", i, x_ant_data_array[i]);
    //    printf("y_ant_data_array[%d] = %f\n", i, y_ant_data_array[i]);
    //    printf("z_ant_data_array[%d] = %f\n", i, z_ant_data_array[i]);
    //    printf("ref_range_data_array[%d] = %f\n", i, ref_range_data_array[i]);
    //}

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

    // Add 1 to accommodate for overlap (helps with interpolation at boundaries)
    //for(int i = 0; i < PULSES*RC_SAMPLES; i++) {
    //    rc_array[i] = (cfloat) {i, i};
    //}
    
    // AZIMUTH RESOLUTION GRID
    double az_res = 0;
    double half_az_width = 0;
    if (PULSES != 1) {
        double az_ant[PULSES];
        for (int i = 0; i < PULSES; i++) {
            az_ant[i] = std::atan2(y_ant_data_array[i], x_ant_data_array[i]);
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

    // TARGET PIXELS
    //cfloat* xy_px_array = (cfloat*) GMIO::malloc(PULSES*RC_SAMPLES*sizeof(cfloat));
    //float* z_px_array = (float*) GMIO::malloc(PULSES*RC_SAMPLES*sizeof(float));
    //for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {
    //    for(int rng_idx = 0; rng_idx < RC_SAMPLES; rng_idx++) {
    //        int idx = pulse_idx*RC_SAMPLES + rng_idx;
    //        xy_px_array[idx] = (cfloat) {(rng_idx-HALF_RANGE_SAMPLES)*RANGE_RES, az_res*pulse_idx - half_az_width};
    //        z_px_array[idx] = 0.0;

    //        printf("pixels[%d] = {%f, %f}\n", idx, xy_px_array[idx].real, xy_px_array[idx].imag);
    //    }
    //}
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

            printf("xyz_px_array[%d] = {%f, %f, %f}\n", idx-3, xyz_px_array[idx-3], xyz_px_array[idx-2], xyz_px_array[idx-1]);
        }
    }
    
    //int img_elem_size = PULSES*RC_SAMPLES;
    int img_byte_size = PULSES*RC_SAMPLES*sizeof(cfloat);
    cfloat* img_array = (cfloat*) GMIO::malloc(img_byte_size);
    //for(int i = 0; i < RC_SAMPLES; i++) {
    //    img_array[i] = (cfloat) {0, 0};
    //}

    //for(int r = 0; r < MAT_ROWS; r++) {
    //    for(int c = 0; c < MAT_COLS; c++) {
    //        int index = (r*MAT_COLS)+c;
    //        float rand_float = static_cast<float>(std::rand()) / RAND_MAX;
    //        gold_data_array[index] = (cfloat) {rand_float, 0};
    //    }
    //}
    
    // Number of iteration for the AIE graphs to run
    for (int inst = 0; inst < INSTANCES; inst++) {
        //fftGraph[inst].run();
        //hpGraph[inst].run();
        //ifftGraph[inst].run();
        ////cplxConjGraph[inst].run();
        bpGraph[inst].run();
    }

    //int fft_per_ssr_entry_size = MAT_COLS / FFT_NPORTS;
    //int fft_per_ssr_byte_size = fft_per_ssr_entry_size * sizeof(cfloat);

    //TODO: DEBUG
    //print_arr(gold_data_array, 2, 5);

    // RTP Params
    //int32 rtp_img_elem_cnt_result[4] = {};
   
    // Loop through pipeline ITER times
    int inst = 0;
    
    for(int iter=0; iter<ITER; iter++) {
        printf("\nPERFORM BACKPROJECTION (ITER = %d) (INST = %d)\n", iter, inst);

        int px_per_ai = (PULSES*RC_SAMPLES)/IMG_SOLVERS;
        int rc_per_ai = RC_SAMPLES/IMG_SOLVERS;
        printf("px_per_ai = %d\n", px_per_ai);
        
        // Pass in slowtime data into AI kernels
        bpGraph[inst].gmio_in_x_ant_pos.gm2aie_nb(x_ant_data_array, PULSES*sizeof(float));
        bpGraph[inst].gmio_in_y_ant_pos.gm2aie_nb(y_ant_data_array, PULSES*sizeof(float));
        bpGraph[inst].gmio_in_z_ant_pos.gm2aie_nb(z_ant_data_array, PULSES*sizeof(float));
        bpGraph[inst].gmio_in_ref_range.gm2aie_nb(ref_range_data_array, PULSES*sizeof(float));

        // Pass in other data into bp AI kernels
        for(int pulse_idx=0; pulse_idx<PULSES; pulse_idx++) {

            bpGraph[inst].gmio_in_rc.gm2aie_nb(rc_array + pulse_idx*RC_SAMPLES, RC_SAMPLES*sizeof(cfloat));
            bpGraph[inst].gmio_in_xyz_px[0].gm2aie_nb(xyz_px_array, PULSES*RC_SAMPLES*sizeof(float)*3);

            for(int kern_id=0; kern_id<IMG_SOLVERS; kern_id++) {
                // Dump image if on last pulse, otherwise keep focusing the image
                if (pulse_idx == PULSES-1) {
                    bpGraph[inst].update(bpGraph[inst].rtp_dump_img_in[kern_id], true);
                } else {
                    bpGraph[inst].update(bpGraph[inst].rtp_dump_img_in[kern_id], false);
                }
                
                // IS IT POSSIBLE TO DO MULTIPLE ITERATIONS OF KERNEL WITHOUT IT COUNTING AS A PULESE, BUT JUST TO FORCE MORE 
                // PIXELS INTO KERNEL FOR CUMULATION? I FEEL THIS FEATURE WOULD BE GOOD FOR EXSTENSABILITY
                //bpGraph[inst].gmio_in_xy_px[kern_id].gm2aie_nb(xy_px_array + kern_id*px_per_ai, px_per_ai*sizeof(cfloat));
                //bpGraph[inst].gmio_in_z_px[kern_id].gm2aie_nb(z_px_array + kern_id*px_per_ai, px_per_ai*sizeof(float));

                //float low_x_px_idx_bound = (range_width_seg*3)-(kern_id*range_width_seg);
                //float high_x_px_idx_bound = (range_width_seg*4)-(kern_id*range_width_seg);

                //float start_range_px = -range_width/2 + range_width_seg*kern_id;
                //float end_range_px = -range_width/2 + range_width_seg*(kern_id+1);
                //float start_az_px = az_width/2.0;
                //float end_az_px = -az_width/2.0;

                //float start_range_px = -range_width/2;
                //float end_range_px = range_width/2;

                //printf("start_range_px=%f | end_range_px=%f | start_az_px=%f | end_az_px=%f | range_res=%f\n", 
                //        start_range_px, end_range_px, start_az_px, end_az_px, range_res);


                //2 PULSES SEEMED VERY CLOSE TO MATLAB. 3 PULSES DID NOT WORK;IT GAVE A WEIRD PIC IN MATLAB. 4 SEEMED OKAY BUT SLIGHTLY DIFFRENT THAN MATLAB. TRYING 10 TO SEE IF THAT WORKS, OR IF IT NEEDS TO BE A FACTOR OF THE POWER OF 2

                // Derive RC offsets for AI kernel per pulse. 
                // The motivation behind the following derivation is to try and pass "usable" range 
                // compressed values into each img_reconstruct_kern AI kernel. Because we calculate the 
                // boundaries of what the AI kernels will use here (on host/ARM), that gives us insight
                // into which range compressed values will actually be utilized for a specific AI kernel 
                // given their target pixels they are deriving for. This reduces the necessary size of the
                // rc_in ping-pong buffer into the kernel (which is necessary because we are already pushing
                // the stack limit). Note: This same calculation occurs on the AI kernel for each target pixel.
                //float upper_x_per_ai = xy_px_array[kern_id*px_per_ai + (px_per_ai-1)].real;
                //float upper_y_per_ai = xy_px_array[kern_id*px_per_ai + (px_per_ai-1)].imag;
                //float upper_z_per_ai = z_px_array[0]; // All Z pixels currently 0, so upper bound is any Z pixel
                //float dR_bounds = sqrt(pow(x_ant_data_array[pulse_idx] - upper_x_per_ai, 2) 
                //                     + pow(y_ant_data_array[pulse_idx] - upper_y_per_ai, 2) 
                //                     + pow(z_ant_data_array[pulse_idx] - upper_z_per_ai, 2))
                //                     - ref_range_data_array[pulse_idx];

                //float px_idx_bound = dR_bounds/range_res + half_range_samples;
                //int rounded_px_idx_bound = (int) std::floor(px_idx_bound);

                //// The offset for gm2aie needs to be 128 bit aligned. Because cfloats are 8B (64b),
                //// then the rc_idx_offset just needs to be even
                //int rc_idx_offset = rounded_px_idx_bound - rounded_px_idx_bound%2;
                //printf("%d: idx: %d | upper_x: %f | upper_y: %f | upper_z: %f | rounded_px_idx_bound: %d | rc_idx_offset: %d\n", 
                //        kern_id, kern_id*px_per_ai + (px_per_ai-1), upper_x_per_ai, upper_y_per_ai, upper_z_per_ai, rounded_px_idx_bound, rc_idx_offset);

                ////printf("dR_bounds: %f | px_idx_bound: %f | rounded_px_idx_bound: %d | rc_idx_offset: %d\n", dR_bounds, px_idx_bound, rounded_px_idx_bound, rc_idx_offset);

                //bpGraph[inst].update(bpGraph[inst].rtp_rc_idx_offset_in[kern_id], rc_idx_offset);
                bpGraph[inst].update(bpGraph[inst].rtp_rc_idx_offset_in[kern_id], 0);


                // Need to be wiser and stratigically pass in data based on what the input target pixels are for that AI tile
                //bpGraph[inst].gmio_in_rc[kern_id].gm2aie_nb(rc_array, RC_SAMPLES*sizeof(cfloat));
                //bpGraph[inst].gmio_in_rc[kern_id].gm2aie_nb(rc_array + (3-kern_id)*px_per_ai, px_per_ai*sizeof(cfloat));
                //bpGraph[inst].gmio_in_rc[kern_id].gm2aie_nb(rc_array + pulse_idx*RC_SAMPLES, RC_SAMPLES*sizeof(cfloat));
                
                bpGraph[inst].gmio_out_img[kern_id].aie2gm_nb(img_array + kern_id*px_per_ai, px_per_ai*sizeof(cfloat));
            }
            for(int kern_id=0; kern_id<IMG_SOLVERS; kern_id++) {
                bpGraph[inst].gmio_out_img[kern_id].wait();
                printf("after wait sig\n");
            }
        }

        //bpGraph[inst].gmio_out_img.aie2gm_nb(img_array, img_byte_size);

        //for(int kern_id=0; kern_id<IMG_SOLVERS; kern_id++) {
        //    // RTP header idx
        //    bpGraph[inst].read(bpGraph[inst].rtp_img_elem_cnt_out[kern_id], rtp_img_elem_cnt_result[kern_id]);
        //    printf("rtp_img_elem_cnt_result[%d] = %d\n", kern_id, rtp_img_elem_cnt_result[kern_id]);
        //}

        //bpGraph[inst].gmio_out_img.wait();
        //printf("after wait sig\n");
        //for(int kern_id=0; kern_id<IMG_SOLVERS; kern_id++) {
        //    bpGraph[inst].gmio_out_img[kern_id].wait();
        //    printf("after wait sig\n");
        //}

        //printf("START TIMER...\n");
        //sleep(120);
        //std::this_thread::sleep_for(std::chrono::seconds(15));
        //printf("DONE TIMER\n");


        //int pkt_header_id = *(unsigned int*)&img_array[0].real & 0x1F;
        //printf("array[%d] = %d\n", 0,  pkt_header_id);

        //TODO: DEBUG
        //printf("array[%d] = 0x%08x\n", 0,  *(unsigned int*)&img_array[0].real);
        //printf("array[%d] = 0x%08x\n", 10, *(unsigned int*)&img_array[10].imag);
        //printf("array[%d] = 0x%08x\n", 30, *(unsigned int*)&img_array[30].real);
        //printf("array[%d] = 0x%08x\n", 31, *(unsigned int*)&img_array[31].real);
        //printf("array[%d] = 0x%08x\n", 32, *(unsigned int*)&img_array[32].real);
        //printf("array[%d] = 0x%08x\n", 56, *(unsigned int*)&img_array[56].imag);


        // Open a file for writing
        // Current working dir inside build/hw/aiesim/
        FILE *img_fp = fopen("../../../design/aie/img.csv", "w");
        if (img_fp == NULL) {
            perror("Error opening img.csv file");
            return EXIT_FAILURE;
        }

        fprintf(img_fp, "%.12f%+.12fi", img_array[0].real, img_array[0].imag);
        for(int i=1; i<PULSES*RC_SAMPLES; i++) {
            if (i%RC_SAMPLES == 0) {
                fprintf(img_fp, "\n");
            }
            fprintf(img_fp, ",%.12f%+.12fi", img_array[i].real, img_array[i].imag);
        }

        for(int i=0; i<PULSES*RC_SAMPLES; i++) {
            printf("array[%d] = {%f, %f}\n", i, img_array[i].real, img_array[i].imag);
        }

        //for(int i=2048; i<2048+32; i++) {
        //    printf("array[%d] = {%f, %f}\n", i, img_array[i].real, img_array[i].imag);
        //}

        //for(int i=4096; i<4096+32; i++) {
        //    printf("array[%d] = {%f, %f}\n", i, img_array[i].real, img_array[i].imag);
        //}

        //for(int i=6144; i<6144+32; i++) {
        //    printf("array[%d] = {%f, %f}\n", i, img_array[i].real, img_array[i].imag);
        //}

        //for(int i=0; i<DIFFER_RANGE_SOLVERS; i++) {
        //    bpGraph[0].gmio_in_xy_px[i].gm2aie_nb(xy_px_array + i*2048, 2048*sizeof(cfloat));
        //}

        //for(int i=0; i<IMG_SOLVERS; i++) {
        //    bpGraph[0].gmio_in_rc[i].gm2aie_nb(rc_array + i*2048, 2048*sizeof(cfloat));
        //    bpGraph[0].gmio_out_img[i].aie2gm_nb(img_array + i*2048, 2048*sizeof(cfloat));
        //}

        //for(int inst=0; inst<INSTANCES; inst++) {
        //    for(int bp=0; bp<IMG_SOLVERS; bp++) {
        //        bpGraph[inst].gmio_out_img[bp].wait();
        //    }
        //}

        //for(int i=0; i<32; i++) {
        //    printf("array[%d] = {%f, %f}\n", i, img_array[i].real, img_array[i].imag);
        //}

        //for(int i=2048; i<2048+32; i++) {
        //    printf("array[%d] = {%f, %f}\n", i, img_array[i].real, img_array[i].imag);
        //}

        //for(int i=4096; i<4096+32; i++) {
        //    printf("array[%d] = {%f, %f}\n", i, img_array[i].real, img_array[i].imag);
        //}

        //for(int i=6144; i<6144+32; i++) {
        //    printf("array[%d] = {%f, %f}\n", i, img_array[i].real, img_array[i].imag);
        //}

        //for(int inst=0; inst<INSTANCES; inst++) {
        //    printf("\nPERFORM FFT (ITER = %d) (INST = %d)\n", iter, inst);

        //    // Set up AIE graph to run MAT_ROWS times
        //    fftGraph[inst].run(MAT_ROWS);

        //    // Pass in data to AIE
        //    for(int row=0; row<MAT_ROWS; row++) {
        //        for(int ssr=0; ssr<FFT_NPORTS; ssr++) {
        //            fftGraph[inst].gmio_in[ssr].gm2aie_nb(gold_data_array + MAT_COLS*row + ssr*fft_per_ssr_entry_size, fft_per_ssr_byte_size);
        //            fftGraph[inst].gmio_out[ssr].aie2gm_nb(data_arrays[inst] + MAT_COLS*row + ssr*fft_per_ssr_entry_size, fft_per_ssr_byte_size);
        //        }
        //    }
        //}

        //// Block until AIE finishes
        //for(int inst=0; inst<INSTANCES; inst++) {
        //    for(int ssr=0; ssr<FFT_NPORTS; ssr++) {
        //        fftGraph[inst].gmio_out[ssr].wait();
        //    }
        //    //TODO: DEBUG
        //    print_arr(data_arrays[inst], MAT_ROWS, MAT_COLS);
        //}

        ////std::cout << "\nPERFORM COMPLEX CONJUGATE" << std::endl;
        ////for(int inst=0; inst<CPLX_CONJ_INSTANCES; inst++) {
        ////    cplxConjGraph[inst].gmio_in.gm2aie_nb(data_array+inst*CPLX_CONJ_POINT_SIZE, BLOCK_SIZE_BYTES/CPLX_CONJ_INSTANCES);
        ////    cplxConjGraph[inst].gmio_out.aie2gm_nb(data_array+inst*CPLX_CONJ_POINT_SIZE, BLOCK_SIZE_BYTES/CPLX_CONJ_INSTANCES);
        ////}

        ////// Block until AIE finishes
        ////for(int inst=0; inst<CPLX_CONJ_INSTANCES; inst++) {
        ////    cplxConjGraph[inst].gmio_out.wait();
        ////}

        //////TODO: DEBUG
        ////print_arr(data_array, 2, MAT_COLS);

        //int hp_per_ssr_byte_size = BLOCK_SIZE_BYTES / TP_SSR;
        //int hp_per_ssr_entry_size = BLOCK_SIZE_ENTRIES / TP_SSR;
        //for(int r = 0; r < MAT_ROWS; r++) {
        //    for(int c = 0; c < MAT_COLS; c++) {
        //        int index = (r*MAT_COLS)+c;
        //        one_array[index] = (cfloat) {1, 0};
        //    }
        //}
        //for(int inst=0; inst<INSTANCES; inst++) {
        //    printf("\nPERFORM ELEMENT-WISE MATRIX MULTIPLY (ITER = %d) (INST = %d)\n", iter, inst);
        //    hpGraph[inst].run(MAT_ROWS/TP_NUM_FRAMES);
        //    for(int ssr=0; ssr<TP_SSR; ssr++) {
        //        hpGraph[inst].gmio_in_A[ssr].gm2aie_nb(one_array+ssr*hp_per_ssr_entry_size, hp_per_ssr_byte_size);
        //        hpGraph[inst].gmio_in_B[ssr].gm2aie_nb(data_arrays[inst]+ssr*hp_per_ssr_entry_size, hp_per_ssr_byte_size);
        //        hpGraph[inst].gmio_out[ssr].aie2gm_nb(data_arrays[inst]+ssr*hp_per_ssr_entry_size, hp_per_ssr_byte_size);
        //    }
        //}

        //// Block until AIE finishes
        //for(int inst=0; inst<INSTANCES; inst++) {
        //    for(int ssr=0; ssr<TP_SSR; ssr++) {
        //        hpGraph[inst].gmio_out[ssr].wait();
        //    }
        //    //TODO: DEBUG
        //    print_arr(data_arrays[inst], MAT_ROWS, MAT_COLS);
        //}

        //for(int inst=0; inst<INSTANCES; inst++) {
        //    printf("\nPERFORM IFFT (ITER = %d) (INST = %d)\n", iter, inst);
        //    ifftGraph[inst].run(MAT_ROWS);
        //    for(int row=0; row<MAT_ROWS; row++) {
        //        for(int ssr=0; ssr<FFT_NPORTS; ssr++) {
        //            ifftGraph[inst].gmio_in[ssr].gm2aie_nb(data_arrays[inst] + MAT_COLS*row + ssr*fft_per_ssr_entry_size, fft_per_ssr_byte_size);
        //            ifftGraph[inst].gmio_out[ssr].aie2gm_nb(data_arrays[inst] + MAT_COLS*row + ssr*fft_per_ssr_entry_size, fft_per_ssr_byte_size);
        //        }
        //    }
        //}

        //// Block until AIE finishes
        //for(int inst=0; inst<INSTANCES; inst++) {
        //    for(int ssr=0; ssr<FFT_NPORTS; ssr++) {
        //        ifftGraph[inst].gmio_out[ssr].wait();
        //    }
        //    total_err_cnt+=ifftErrorCheck(data_arrays[inst], gold_data_array);
        //}
    }

    // End AIE graphs
    //for(int inst=0; inst<INSTANCES; inst++) {
    //    fftGraph[inst].end();
    //    ifftGraph[inst].end();
    //    hpGraph[inst].end();
    //}
    //for(int inst=0; inst<CPLX_CONJ_INSTANCES*INSTANCES; inst++) {
    //    cplxConjGraph[inst].end();
    //}

    return total_err_cnt;
}

//const int BLOCK_SIZE_ENTRIES = MAT_ROWS * MAT_COLS;
//const int BLOCK_SIZE_BYTES = BLOCK_SIZE_ENTRIES * sizeof(cfloat);
//
//int fftRowErrorCheck(cfloat* fft_row_array, int instance = -1) {
//    int error_cnt = 0;
//
//    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
//        if(i < MAT_ROWS) {
//            if(fft_row_array[i].real != FFT_SAMPLE_DATA || fft_row_array[i].imag != FFT_SAMPLE_DATA)
//                error_cnt++;
//        } else {
//            if(fft_row_array[i].real != 0 || fft_row_array[i].imag != 0)
//                error_cnt++;
//        }
//    }
//
//    if (instance >= 0)
//        printf("INSTANCE %d: ", instance);
//    printf("FFT ROW ERRORS: %d\n", error_cnt);
//
//    return error_cnt;
//}
//
//int fftColErrorCheck(cfloat* fft_col_array, int instance = -1) {
//    int error_cnt = 0;
//
//    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
//        if(fft_col_array[i].real != FFT_SAMPLE_DATA || fft_col_array[i].imag != FFT_SAMPLE_DATA)
//            error_cnt++;
//    }
//
//    if (instance >= 0)
//        printf("INSTANCE %d: ", instance);
//    printf("FFT COL ERRORS: %d\n", error_cnt);
//
//    return error_cnt;
//}
//
//int cplxConjErrorCheck(cfloat* cplx_conj_array, int instance = -1) {
//    int error_cnt = 0;
//
//    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
//        if(cplx_conj_array[i].real != FFT_SAMPLE_DATA || cplx_conj_array[i].imag != -FFT_SAMPLE_DATA)
//            error_cnt++;
//    }
//
//    if (instance >= 0)
//        printf("INSTANCE %d: ", instance);
//    printf("COMPLEX CONJUGATE ERRORS: %d\n", error_cnt);
//
//    return error_cnt;
//}
//
//int hpErrorCheck(cfloat* hp_array, int instance = -1) {
//    int error_cnt = 0;
//
//    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
//        if(hp_array[i].real != HP_SAMPLE_DATA || hp_array[i].imag != 0)
//            error_cnt++;
//    }
//
//    if (instance >= 0)
//        printf("INSTANCE %d: ", instance);
//    printf("ELEMENT-WISE MATRIX MULTIPLY ERRORS: %d\n", error_cnt);
//
//    return error_cnt;
//}
//
//int ifftColErrorCheck(cfloat* ifft_col_array, int instance = -1) {
//    int error_cnt = 0;
//
//    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
//        // First column should be HP_SAMPLE_DATA*MAT_COLS and everywhere else 0
//        if(i%MAT_COLS == 0) {
//            if(ifft_col_array[i].real != HP_SAMPLE_DATA*MAT_COLS || ifft_col_array[i].imag != 0)
//                error_cnt++;
//        } else {
//            if(ifft_col_array[i].real != 0 || ifft_col_array[i].imag != 0)
//                error_cnt++;
//        }
//    }
//
//    if (instance >= 0)
//        printf("INSTANCE %d: ", instance);
//    printf("IFFT COL ERRORS: %d\n", error_cnt);
//
//    return error_cnt;
//}
//
//int ifftRowErrorCheck(cfloat* ifft_row_array, int instance = -1) {
//    int error_cnt = 0;
//
//    if(ifft_row_array[0].real != HP_SAMPLE_DATA*BLOCK_SIZE_ENTRIES || ifft_row_array[0].imag != 0)
//        error_cnt++;
//
//    for(int i = 1; i < BLOCK_SIZE_ENTRIES; i++) {
//        if(ifft_row_array[i].real != 0 || ifft_row_array[i].imag != 0)
//            error_cnt++;
//    }
//
//    if (instance >= 0)
//        printf("INSTANCE %d: ", instance);
//    printf("IFFT ROW ERRORS: %d\n", error_cnt);
//
//    return error_cnt;
//}
//
//int peakSearchErrorCheck(cfloat* peak_array, int instance = -1) {
//    int error_cnt = 0;
//
//    int max_val = peak_array[0].real;
//    int index = peak_array[0].imag;
//
//    if(max_val != HP_SAMPLE_DATA*BLOCK_SIZE_ENTRIES || index != 0)
//        error_cnt++;
//
//    if (instance >= 0)
//        printf("INSTANCE %d: ", instance);
//    printf("PEAK SEARCH ERRORS: %d\n", error_cnt);
//
//    return error_cnt;
//}
//
//
//int main(int argc, char ** argv) {
//
//    // Iterations to run. If you want to increment this, make sure there is enough sim data 
//    // within the ./aie/aiesim_data/ directory to simulate multiple iteratings. You can do
//    // this by matching this number with the ITER_CNT number in the sq_matrix_input_generate.py
//    // script and running that script to generate the proper amount of data.
//    const int ITER = 1;
//
//    // Total error count
//    int total_err_cnt = 0;
//
//    // Initialize all AIE graphs
//    for(int inst=0; inst<INSTANCES; inst++) {
//        fftRowsGraph[inst].init();
//        fftColsGraph[inst].init();
//        ifftColsGraph[inst].init();
//        ifftRowsGraph[inst].init();
//        cplxConjGraph[inst].init();
//        hpGraph[inst].init();
//        peakGraph[inst].init();
//    }
//    
//    // Create a map_fft_array placeholder so we don't need to keep computing the 
//    // map fft every time
//    cfloat* map_fft_array = (cfloat*) GMIO::malloc(BLOCK_SIZE_BYTES);
//    
//    // Create and populate the map image with input data that can be validated
//    // throughout the pipelines
//    cfloat* map_array = (cfloat*) GMIO::malloc(BLOCK_SIZE_BYTES);
//    map_array[0] = (cfloat) {FFT_SAMPLE_DATA, FFT_SAMPLE_DATA};
//    for(int j = 1; j < BLOCK_SIZE_ENTRIES; j++) {
//        map_array[j] = (cfloat) {0, 0};
//    }
//    
//    // Malloc space for template images
//    std::vector<cfloat*> tmpl_arrays;
//    for(int inst=0; inst<INSTANCES; inst++) {
//        tmpl_arrays.push_back((cfloat*) GMIO::malloc(BLOCK_SIZE_BYTES));
//    }
//    
//    // Number of iteration for the AIE graphs to run
//    fftRowsGraph[0].run(MAT_ROWS + (MAT_ROWS*ITER));
//    fftColsGraph[0].run(MAT_COLS + (MAT_COLS*ITER));
//    ifftColsGraph[0].run(MAT_COLS*ITER);
//    ifftRowsGraph[0].run(MAT_ROWS*ITER);
//    cplxConjGraph[0].run(MAT_ROWS*ITER);
//    hpGraph[0].run((TP_DIM/TP_NUM_FRAMES)*ITER);
//    peakGraph[0].run(ITER);
//    for(int inst=1; inst<INSTANCES; inst++) {
//        fftRowsGraph[inst].run(MAT_ROWS*ITER);
//        fftColsGraph[inst].run(MAT_COLS*ITER);
//        ifftColsGraph[inst].run(MAT_COLS*ITER);
//        ifftRowsGraph[inst].run(MAT_ROWS*ITER);
//        cplxConjGraph[inst].run(MAT_ROWS*ITER);
//        hpGraph[inst].run((TP_DIM/TP_NUM_FRAMES)*ITER);
//        peakGraph[inst].run(ITER);
//    }
//
//
//    std::cout << "\nPERFORM 2D FFT OF MAP IMAGE" << std::endl;
//    // Perform row-wise 1D FFT on map image
//    fftRowsGraph[0].gmio_in.gm2aie_nb(map_array, BLOCK_SIZE_BYTES);
//    fftRowsGraph[0].gmio_out.aie2gm_nb(map_fft_array, BLOCK_SIZE_BYTES);
//    fftRowsGraph[0].gmio_out.wait();
//    total_err_cnt+=fftRowErrorCheck(map_fft_array);
//
//    // Perform col-wise 1D FFT of map image (PL sends AIE the transposed data)
//    fftColsGraph[0].gmio_out.aie2gm_nb(map_fft_array, BLOCK_SIZE_BYTES);
//    fftColsGraph[0].gmio_out.wait();
//    total_err_cnt+=fftColErrorCheck(map_fft_array);
//
//    // Start image cross correlations between the map and template images for ITER iterations
//    for(int n=0; n<ITER; n++) {
//        // Create and populate the template images with input data that can be validated
//        // throughout the pipeline
//        for(int inst=0; inst<INSTANCES; inst++) {
//            tmpl_arrays[inst][0] = (cfloat) {FFT_SAMPLE_DATA, FFT_SAMPLE_DATA};
//            for(int j = 1; j < BLOCK_SIZE_ENTRIES; j++) {
//                tmpl_arrays[inst][j] = (cfloat) {0, 0};
//            }
//        }
//
//
//        std::cout << "\nPERFORM 2D FFT OF TEMPLATE IMAGES" << std::endl;
//        // Perform row-wise 1D FFT on all tmplate images at the same time
//        for(int inst=0; inst<INSTANCES; inst++) {
//            fftRowsGraph[inst].gmio_in.gm2aie_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
//            fftRowsGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
//        }
//
//        // Block until AIE sends host computed 1D FFT row-wise results
//        for(int inst=0; inst<INSTANCES; inst++) {
//            fftRowsGraph[inst].gmio_out.wait();
//            total_err_cnt+=fftRowErrorCheck(tmpl_arrays[inst], inst);
//        }
//
//        for(int inst=0; inst<INSTANCES; inst++) {
//        }
//
//        // Perform col-wise 1D FFT on all template images (PL sends AIE the transposed data)
//        for(int inst=0; inst<INSTANCES; inst++) {
//            fftColsGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
//        }
//
//        // Block until AIE finishes
//        for(int inst=0; inst<INSTANCES; inst++) {
//            fftColsGraph[inst].gmio_out.wait();
//            total_err_cnt+=fftColErrorCheck(tmpl_arrays[inst], inst);
//        }
//
//
//        std::cout << "\nPERFORM COMPLEX CONJUGATE ON TEMPLATE IMAGES" << std::endl;
//        // Perform complex conjugate on all tmplate images at the same time
//        for(int inst=0; inst<INSTANCES; inst++) {
//            cplxConjGraph[inst].gmio_in.gm2aie_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
//            cplxConjGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
//        }
//
//        // Block until AIE finishes
//        for(int inst=0; inst<INSTANCES; inst++) {
//            cplxConjGraph[inst].gmio_out.wait();
//            total_err_cnt+=cplxConjErrorCheck(tmpl_arrays[inst], inst);
//        }
//
//        std::cout << "\nPERFORM ELEMENT-WISE MATRIX MULTIPLY BETWEEN MAP AND TEMPLATE IMAGES" << std::endl;
//        // Perform element-wise matrix multiply te on all tmplate images at the same time
//        int per_ssr_byte_size = BLOCK_SIZE_BYTES / TP_SSR;
//        int per_ssr_entry_size = BLOCK_SIZE_ENTRIES / TP_SSR;
//        for(int inst=0; inst<INSTANCES; inst++) {
//            for(int ssr=0; ssr<TP_SSR; ssr++) {
//                hpGraph[inst].gmio_in_A[ssr].gm2aie_nb(map_fft_array, per_ssr_byte_size);
//                hpGraph[inst].gmio_in_B[ssr].gm2aie_nb(tmpl_arrays[inst], per_ssr_byte_size);
//                hpGraph[inst].gmio_out[ssr].aie2gm_nb(tmpl_arrays[inst]+ssr*per_ssr_entry_size, per_ssr_byte_size);
//            }
//        }
//
//        // Block until AIE finishes
//        for(int inst=0; inst<INSTANCES; inst++) {
//            for(int ssr=0; ssr<TP_SSR; ssr++) {
//                hpGraph[inst].gmio_out[ssr].wait();
//            }
//            total_err_cnt+=hpErrorCheck(tmpl_arrays[inst], inst);
//        }
//
//
//        std::cout << "\nPERFORM 2D iFFT OF IMAGES" << std::endl;
//        // Perform col-wise 1D iFFT on all images at the same time
//        for(int inst=0; inst<INSTANCES; inst++) {
//            ifftColsGraph[inst].gmio_in.gm2aie_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
//            ifftColsGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
//        }
//
//        // Block until AIE finishes
//        for(int inst=0; inst<INSTANCES; inst++) {
//            ifftColsGraph[inst].gmio_out.wait();
//            total_err_cnt+=ifftColErrorCheck(tmpl_arrays[inst], inst);
//        }
//
//        // Perform row-wise 1D iFFT on all images (PL sends AIE the transposed data)
//        for(int inst=0; inst<INSTANCES; inst++) {
//            ifftRowsGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
//        }
//
//        // Block until AIE finishes
//        for(int inst=0; inst<INSTANCES; inst++) {
//            ifftRowsGraph[inst].gmio_out.wait();
//            total_err_cnt+=ifftRowErrorCheck(tmpl_arrays[inst], inst);
//        }
//
//
//        std::cout << "\nPERFORM PEAK SEARCH ON ALL IMAGES" << std::endl;
//        // Perform peak search on all images at the same time
//        for(int inst=0; inst<INSTANCES; inst++) {
//            peakGraph[inst].gmio_in.gm2aie_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
//            peakGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], sizeof(cfloat));
//        }
//
//        // Block until AIE finishes
//        for(int inst=0; inst<INSTANCES; inst++) {
//            peakGraph[inst].gmio_out.wait();
//            total_err_cnt+=peakSearchErrorCheck(tmpl_arrays[inst], inst);
//        }
//
//    }
//
//    if(total_err_cnt) {
//        printf("\nTEST FAILED!!!\n\n");
//    } else {
//        printf("\nTEST PASSED!!!\n\n");
//    }
//
//    // Free memory
//    GMIO::free(map_array);
//    GMIO::free(map_fft_array);
//    for(int inst=0; inst<INSTANCES; inst++) {
//        GMIO::free(tmpl_arrays[inst]);
//
//        // End AIE graphs
//        fftRowsGraph[inst].end();
//        fftColsGraph[inst].end();
//        ifftRowsGraph[inst].end();
//        ifftColsGraph[inst].end();
//        cplxConjGraph[inst].end();
//        hpGraph[inst].end();
//        peakGraph[inst].end();
//    }
//
//    return total_err_cnt;
//}

#endif
