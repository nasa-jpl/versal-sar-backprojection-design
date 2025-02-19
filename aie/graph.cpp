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

const float EPSILON = 0.0001;
const int BLOCK_SIZE_ENTRIES = MAT_ROWS*MAT_COLS;
const int BLOCK_SIZE_BYTES = BLOCK_SIZE_ENTRIES * sizeof(TT_DATA);

// FFT output should all be FFT_SAMPLE_DATA for both real and imaginary
//int fftErrorCheck(TT_DATA* data_array) {
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

int ifftErrorCheck(TT_DATA* data_array, TT_DATA* gold_data_array) {
    int error_cnt = 0;
    float real, imag;
    float gold_real, gold_imag;

    for(int r = 0; r < MAT_ROWS; r++) {
        for(int c = 0; c < MAT_COLS; c++) {
            int index = (r*MAT_COLS)+c;
            real = data_array[index].real/TP_POINT_SIZE;
            imag = data_array[index].imag/TP_POINT_SIZE;
            gold_real = gold_data_array[index].real;
            //gold_imag = -gold_data_array[index].imag; // Adding negative to accommodate for complx conjugation
            gold_imag = gold_data_array[index].imag; // Adding negative to accommodate for complx conjugation
            //printf("IFFT[%d, %d] = {%f, %f}\n", r, c, real, imag);

            if(std::abs(real-gold_real) > EPSILON || std::abs(imag-gold_imag) > EPSILON) {
                printf("ERROR IFFT[%d, %d] = {%f, %f} | EXPECTED {%f, %f}\n", r, c, real, imag, gold_real, gold_imag);
                error_cnt++;
            }
        }
    }

    printf("IFFT ERRORS: %d\n\n", error_cnt);

    return error_cnt;
}

void print_arr(TT_DATA* arr, int num_rows, int num_cols) {
    for(int r=0; r<num_rows; r++) {
        for(int c=0; c<num_cols; c++) {
            int index = (r*MAT_COLS) + c;
            printf("array[%d][%d] = {%f, %f}\n", r, c, arr[index].real, arr[index].imag);
        }
    }
}

void reorderDataArray(TT_DATA* data_array) {
    // Create a temporary array to hold the rearranged data
    TT_DATA* temp = (TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES);
    
    // Rearrange the data into the temp array
    int chunk_size = MAT_COLS;
    int num_datasets = BLOCK_SIZE_ENTRIES / chunk_size;

    for (int dataset = 0; dataset < num_datasets; dataset++) {
        int offset = dataset * chunk_size;

        for (int i = 0; i < chunk_size; i++) {
            int new_index = offset + (i / 4) + (i % 4) * (chunk_size / 4);
            temp[new_index] = data_array[offset + i];
        }
    }

    // Copy the rearranged data back into the original array
    for (int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
        data_array[i] = temp[i];
    }

    printf("\nFLIPPED DATA_ARRAY\n");
    for(int r=0; r<MAT_ROWS; r++) {
        for(int c=0; c<MAT_COLS; c++) {
            int index = (r*MAT_COLS) + c;
            printf("data_array[%d][%d] = {%f, %f}\n", r, c, data_array[index].real, data_array[index].imag);
        }
    }

    // Free the temporary array
    free(temp);
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
    //std::vector<TT_DATA*> data_arrays;
    //for(int inst=0; inst<INSTANCES; inst++) {
    //    data_arrays.push_back((TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES));
    //}
    //TT_DATA* one_array = (TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES);

    // Allocate and populate memory
    //TT_DATA* gold_data_array = (TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES);
    //for(int r = 0; r < MAT_ROWS; r++) {
    //    for(int c = 0; c < MAT_COLS; c++) {
    //        int index = (r*MAT_COLS)+c;
    //        if (c == 1) {
    //            gold_data_array[index] = (TT_DATA) {FFT_SAMPLE_DATA, FFT_SAMPLE_DATA};
    //        }
    //        else {
    //            gold_data_array[index] = (TT_DATA) {0, 0};
    //        }
    //    }
    //}

    // OPEN SAR DATASET FILES
    int PULSES = 2;

    float* x_ant_data_array = (float*) GMIO::malloc(PULSES*sizeof(float));
    float* y_ant_data_array = (float*) GMIO::malloc(PULSES*sizeof(float));
    float* z_ant_data_array = (float*) GMIO::malloc(PULSES*sizeof(float));
    float* ref_range_data_array = (float*) GMIO::malloc(PULSES*sizeof(float));
    
    // File referenced from inside build/hw/aiesim/
    std::ifstream st_file("../../../design/aie/test_data/gotcha_slowtime_pass1_360deg_HH.csv"); // Referenced from build/hw/aiesim/
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

    TT_DATA* rc_array = (TT_DATA*) GMIO::malloc(PULSES*TP_POINT_SIZE*sizeof(TT_DATA));
     
    // File referenced from inside build/hw/aiesim/
    std::string rc_filename = "../../../design/aie/test_data/gotcha_phdata_" + 
                              std::to_string(TP_POINT_SIZE) + 
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
        while (std::getline(ss, value, ',') && rc_samp_cnt < TP_POINT_SIZE) {
            std::regex complex_regex(R"(([+-]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)([+-](?:\d+(?:\.\d*)?|\.\d+)(?:[eE][+-]?\d+)?)i)");
            std::smatch match;
            std::regex_search(value, match, complex_regex);

            float real_part = std::stof(match[1].str());
            float imag_part = std::stof(match[2].str());

            // Store the complex number in the array
            rc_array[pulse_idx*TP_POINT_SIZE + rc_samp_cnt] = (TT_DATA) {real_part, imag_part};

            rc_samp_cnt++;
        }
        rc_samp_cnt = 0;
        pulse_idx++;
    }

    // Add 1 to accommodate for overlap (helps with interpolation at boundaries)
    //for(int i = 0; i < PULSES*TP_POINT_SIZE; i++) {
    //    rc_array[i] = (TT_DATA) {i, i};
    //}

    TT_DATA* xy_px_array = (TT_DATA*) GMIO::malloc(TP_POINT_SIZE*sizeof(TT_DATA));
    float* z_px_array = (float*) GMIO::malloc(TP_POINT_SIZE*sizeof(float));
    const float C = 299792458.0;
    float range_freq_step = 1471301.6;
    //float range_freq_step = 5000000.0;
    int half_range_samples = TP_POINT_SIZE/2;
    float min_freq = 9288080400.0;

    float range_width = C/(2.0*range_freq_step);
    float range_width_seg = range_width/IMG_SOLVERS;
    float range_res = range_width/TP_POINT_SIZE;

    //float total_az = 0.01726837;
    //float az_res = C/(2.0*total_az*min_freq);
    //float az_width = AZ_POINT_SIZE/az_res;

    for(int i = 0; i < TP_POINT_SIZE; i++) {
        xy_px_array[i] = (TT_DATA) {(i-half_range_samples)*range_res, 0};
        z_px_array[i] = 0.0;
    }
    
    int img_elem_size = AZ_POINT_SIZE*TP_POINT_SIZE;
    int img_byte_size = img_elem_size*sizeof(TT_DATA);
    TT_DATA* img_array = (TT_DATA*) GMIO::malloc(img_byte_size);
    //for(int i = 0; i < TP_POINT_SIZE; i++) {
    //    img_array[i] = (TT_DATA) {0, 0};
    //}

    //for(int r = 0; r < MAT_ROWS; r++) {
    //    for(int c = 0; c < MAT_COLS; c++) {
    //        int index = (r*MAT_COLS)+c;
    //        float rand_float = static_cast<float>(std::rand()) / RAND_MAX;
    //        gold_data_array[index] = (TT_DATA) {rand_float, 0};
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

    int fft_per_ssr_entry_size = MAT_COLS / FFT_NPORTS;
    int fft_per_ssr_byte_size = fft_per_ssr_entry_size * sizeof(TT_DATA);

    //TODO: DEBUG
    //print_arr(gold_data_array, 2, 5);

    // RTP Params
    //int32 rtp_img_elem_cnt_result[4] = {};
   
    // Loop through pipeline ITER times
    int inst = 0;
    
    for(int iter=0; iter<ITER; iter++) {
        printf("\nPERFORM BACKPROJECTION (ITER = %d) (INST = %d)\n", iter, inst);

        int tp_pt_sz_per_ai = TP_POINT_SIZE/IMG_SOLVERS;
        printf("tp_pt_sz_per_ai = %d\n", tp_pt_sz_per_ai);
        
        // Pass in slowtime data into AI kernels
        bpGraph[inst].gmio_in_x_ant_pos.gm2aie_nb(x_ant_data_array, PULSES*sizeof(float));
        bpGraph[inst].gmio_in_y_ant_pos.gm2aie_nb(y_ant_data_array, PULSES*sizeof(float));
        bpGraph[inst].gmio_in_z_ant_pos.gm2aie_nb(z_ant_data_array, PULSES*sizeof(float));
        bpGraph[inst].gmio_in_ref_range.gm2aie_nb(ref_range_data_array, PULSES*sizeof(float));

        // Pass in other data into bp AI kernels
        for(int pulse_idx=0; pulse_idx<PULSES; pulse_idx++) {
            for(int kern_id=0; kern_id<IMG_SOLVERS; kern_id++) {
                // Dump image if on last pulse, otherwise keep focusing the image
                if (pulse_idx == PULSES-1) {
                    bpGraph[inst].update(bpGraph[inst].rtp_dump_img_in[kern_id], true);
                } else {
                    bpGraph[inst].update(bpGraph[inst].rtp_dump_img_in[kern_id], false);
                }

                bpGraph[inst].gmio_in_xy_px[kern_id].gm2aie_nb(xy_px_array + kern_id*tp_pt_sz_per_ai, tp_pt_sz_per_ai*sizeof(TT_DATA));
                bpGraph[inst].gmio_in_z_px[kern_id].gm2aie_nb(z_px_array + kern_id*tp_pt_sz_per_ai, tp_pt_sz_per_ai*sizeof(float));
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

                // Derive RC offsets for AI kernel per pulse. 
                // The motivation behind the following derivation is to try and pass "usable" range 
                // compressed values into each img_reconstruct_kern AI kernel. Because we calculate the 
                // boundaries of what the AI kernels will use here (on host/ARM), that gives us insight
                // into which range compressed values will actually be utilized for a specific AI kernel 
                // given their target pixels they are deriving for. This reduces the necessary size of the
                // rc_in ping-pong buffer into the kernel (which is necessary because we are already pushing
                // the stack limit). Note: This same calculation occurs on the AI kernel for each target pixel.
                float dR_bounds = sqrt(pow(x_ant_data_array[pulse_idx]-xy_px_array[kern_id*(tp_pt_sz_per_ai) + (tp_pt_sz_per_ai-1)].real, 2) 
                        + pow(y_ant_data_array[pulse_idx]-xy_px_array[kern_id*(tp_pt_sz_per_ai) + (tp_pt_sz_per_ai-1)].imag, 2) 
                        + pow(z_ant_data_array[pulse_idx]-z_px_array[kern_id*(tp_pt_sz_per_ai) + (tp_pt_sz_per_ai-1)], 2))
                        - ref_range_data_array[pulse_idx];

                float px_idx_bound = dR_bounds/range_res + half_range_samples;
                int rounded_px_idx_bound = (int) std::floor(px_idx_bound);

                // The offset for gm2aie needs to be 128 bit aligned. Because cfloats are 8B (64b),
                // then the rc_idx_offset just needs to be even
                int rc_idx_offset = rounded_px_idx_bound - rounded_px_idx_bound%2;

                bpGraph[inst].update(bpGraph[inst].rtp_rc_idx_offset_in[kern_id], rc_idx_offset);


                // Need to be wiser and stratigically pass in data based on what the input target pixels are for that AI tile
                //bpGraph[inst].gmio_in_rc[kern_id].gm2aie_nb(rc_array, TP_POINT_SIZE*sizeof(TT_DATA));
                //bpGraph[inst].gmio_in_rc[kern_id].gm2aie_nb(rc_array + (3-kern_id)*tp_pt_sz_per_ai, tp_pt_sz_per_ai*sizeof(TT_DATA));
                bpGraph[inst].gmio_in_rc[kern_id].gm2aie_nb(rc_array + pulse_idx*TP_POINT_SIZE + rc_idx_offset, tp_pt_sz_per_ai*sizeof(TT_DATA));
                
                bpGraph[inst].gmio_out_img[kern_id].aie2gm_nb(img_array + kern_id*tp_pt_sz_per_ai, tp_pt_sz_per_ai*sizeof(TT_DATA));
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
        for(int i=0; i<img_elem_size; i++) {
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
        //    bpGraph[0].gmio_in_xy_px[i].gm2aie_nb(xy_px_array + i*2048, 2048*sizeof(TT_DATA));
        //}

        //for(int i=0; i<IMG_SOLVERS; i++) {
        //    bpGraph[0].gmio_in_rc[i].gm2aie_nb(rc_array + i*2048, 2048*sizeof(TT_DATA));
        //    bpGraph[0].gmio_out_img[i].aie2gm_nb(img_array + i*2048, 2048*sizeof(TT_DATA));
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
        //        one_array[index] = (TT_DATA) {1, 0};
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
//const int BLOCK_SIZE_BYTES = BLOCK_SIZE_ENTRIES * sizeof(TT_DATA);
//
//int fftRowErrorCheck(TT_DATA* fft_row_array, int instance = -1) {
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
//int fftColErrorCheck(TT_DATA* fft_col_array, int instance = -1) {
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
//int cplxConjErrorCheck(TT_DATA* cplx_conj_array, int instance = -1) {
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
//int hpErrorCheck(TT_DATA* hp_array, int instance = -1) {
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
//int ifftColErrorCheck(TT_DATA* ifft_col_array, int instance = -1) {
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
//int ifftRowErrorCheck(TT_DATA* ifft_row_array, int instance = -1) {
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
//int peakSearchErrorCheck(TT_DATA* peak_array, int instance = -1) {
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
//    TT_DATA* map_fft_array = (TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES);
//    
//    // Create and populate the map image with input data that can be validated
//    // throughout the pipelines
//    TT_DATA* map_array = (TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES);
//    map_array[0] = (TT_DATA) {FFT_SAMPLE_DATA, FFT_SAMPLE_DATA};
//    for(int j = 1; j < BLOCK_SIZE_ENTRIES; j++) {
//        map_array[j] = (TT_DATA) {0, 0};
//    }
//    
//    // Malloc space for template images
//    std::vector<TT_DATA*> tmpl_arrays;
//    for(int inst=0; inst<INSTANCES; inst++) {
//        tmpl_arrays.push_back((TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES));
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
//            tmpl_arrays[inst][0] = (TT_DATA) {FFT_SAMPLE_DATA, FFT_SAMPLE_DATA};
//            for(int j = 1; j < BLOCK_SIZE_ENTRIES; j++) {
//                tmpl_arrays[inst][j] = (TT_DATA) {0, 0};
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
//            peakGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], sizeof(TT_DATA));
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
