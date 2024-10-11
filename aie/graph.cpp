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
//uint8_t cplx_conj_graph_insts = 0;
uint8_t hp_graph_insts = 0;

const int INSTANCES = 1;
FFTGraph<0,0> fftGraph[INSTANCES];
IFFTGraph<0,2> ifftGraph[INSTANCES];
//CplxConjGraph cplxConjGraph[CPLX_CONJ_INSTANCES*INSTANCES];
HPGraph hpGraph[INSTANCES];

#if defined(__AIESIM__) || defined(__X86SIM__)

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
    const int ITER = 2;

    // Generate random seed
    std::time_t seed = std::time(0);
    std::srand(seed);
    std::cout << "Rand Seed: " << seed << std::endl;

    // Total error count
    int total_err_cnt = 0;

    // Initialize all AIE graphs
    for(int inst=0; inst<INSTANCES; inst++) {
        fftGraph[inst].init();
        hpGraph[inst].init();
        ifftGraph[inst].init();
    }
    //for(int inst=0; inst<CPLX_CONJ_INSTANCES*INSTANCES; inst++) {
    //    cplxConjGraph[inst].init();
    //}
   
    // Allocating memory
    TT_DATA* gold_data_array = (TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES);
    std::vector<TT_DATA*> data_arrays;
    for(int inst=0; inst<INSTANCES; inst++) {
        data_arrays.push_back((TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES));
    }
    TT_DATA* one_array = (TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES);

    // Populating inital sample data
    for(int r = 0; r < MAT_ROWS; r++) {
        for(int c = 0; c < MAT_COLS; c++) {
            int index = (r*MAT_COLS)+c;
            if (c == 1) {
                gold_data_array[index] = (TT_DATA) {FFT_SAMPLE_DATA, FFT_SAMPLE_DATA};
            }
            else {
                gold_data_array[index] = (TT_DATA) {0, 0};
            }
        }
    }

    //for(int r = 0; r < MAT_ROWS; r++) {
    //    for(int c = 0; c < MAT_COLS; c++) {
    //        int index = (r*MAT_COLS)+c;
    //        float rand_float = static_cast<float>(std::rand()) / RAND_MAX;
    //        gold_data_array[index] = (TT_DATA) {rand_float, 0};
    //    }
    //}
    
    // Number of iteration for the AIE graphs to run
    for (int inst = 0; inst < INSTANCES; inst++) {
        fftGraph[inst].run();
        hpGraph[inst].run();
        ifftGraph[inst].run();
    }
    //for (int inst = 0; inst < INSTANCES; inst++) {
    //    fftGraph[inst].run(MAT_ROWS*ITER);
    //    hpGraph[inst].run((MAT_ROWS/TP_NUM_FRAMES)*ITER);
    //    ifftGraph[inst].run(MAT_ROWS*ITER);
    //}

    //for (int inst = 0; inst < CPLX_CONJ_INSTANCES*INSTANCES; inst++) {
    //    cplxConjGraph[inst].run(MAT_ROWS*ITER);
    //}

    int fft_per_ssr_entry_size = MAT_COLS / FFT_NPORTS;
    int fft_per_ssr_byte_size = fft_per_ssr_entry_size * sizeof(TT_DATA);

    //TODO: DEBUG
    print_arr(gold_data_array, 2, 5);
   
    // Loop through pipeline ITER times
    for(int iter=0; iter<ITER; iter++) {
        for(int inst=0; inst<INSTANCES; inst++) {
            printf("\nPERFORM FFT (ITER = %d) (INST = %d)\n", iter, inst);

            // Set up AIE graph to run MAT_ROWS times
            fftGraph[inst].run(MAT_ROWS);

            // Pass in data to AIE
            for(int row=0; row<MAT_ROWS; row++) {
                for(int ssr=0; ssr<FFT_NPORTS; ssr++) {
                    fftGraph[inst].gmio_in[ssr].gm2aie_nb(gold_data_array + MAT_COLS*row + ssr*fft_per_ssr_entry_size, fft_per_ssr_byte_size);
                    fftGraph[inst].gmio_out[ssr].aie2gm_nb(data_arrays[inst] + MAT_COLS*row + ssr*fft_per_ssr_entry_size, fft_per_ssr_byte_size);
                }
            }
        }

        // Block until AIE finishes
        for(int inst=0; inst<INSTANCES; inst++) {
            for(int ssr=0; ssr<FFT_NPORTS; ssr++) {
                fftGraph[inst].gmio_out[ssr].wait();
            }
            //TODO: DEBUG
            print_arr(data_arrays[inst], MAT_ROWS, MAT_COLS);
        }

        //std::cout << "\nPERFORM COMPLEX CONJUGATE" << std::endl;
        //for(int inst=0; inst<CPLX_CONJ_INSTANCES; inst++) {
        //    cplxConjGraph[inst].gmio_in.gm2aie_nb(data_array+inst*CPLX_CONJ_POINT_SIZE, BLOCK_SIZE_BYTES/CPLX_CONJ_INSTANCES);
        //    cplxConjGraph[inst].gmio_out.aie2gm_nb(data_array+inst*CPLX_CONJ_POINT_SIZE, BLOCK_SIZE_BYTES/CPLX_CONJ_INSTANCES);
        //}

        //// Block until AIE finishes
        //for(int inst=0; inst<CPLX_CONJ_INSTANCES; inst++) {
        //    cplxConjGraph[inst].gmio_out.wait();
        //}

        ////TODO: DEBUG
        //print_arr(data_array, 2, MAT_COLS);

        int hp_per_ssr_byte_size = BLOCK_SIZE_BYTES / TP_SSR;
        int hp_per_ssr_entry_size = BLOCK_SIZE_ENTRIES / TP_SSR;
        for(int r = 0; r < MAT_ROWS; r++) {
            for(int c = 0; c < MAT_COLS; c++) {
                int index = (r*MAT_COLS)+c;
                one_array[index] = (TT_DATA) {1, 0};
            }
        }
        for(int inst=0; inst<INSTANCES; inst++) {
            printf("\nPERFORM ELEMENT-WISE MATRIX MULTIPLY (ITER = %d) (INST = %d)\n", iter, inst);
            hpGraph[inst].run(MAT_ROWS/TP_NUM_FRAMES);
            for(int ssr=0; ssr<TP_SSR; ssr++) {
                hpGraph[inst].gmio_in_A[ssr].gm2aie_nb(one_array+ssr*hp_per_ssr_entry_size, hp_per_ssr_byte_size);
                hpGraph[inst].gmio_in_B[ssr].gm2aie_nb(data_arrays[inst]+ssr*hp_per_ssr_entry_size, hp_per_ssr_byte_size);
                hpGraph[inst].gmio_out[ssr].aie2gm_nb(data_arrays[inst]+ssr*hp_per_ssr_entry_size, hp_per_ssr_byte_size);
            }
        }

        // Block until AIE finishes
        for(int inst=0; inst<INSTANCES; inst++) {
            for(int ssr=0; ssr<TP_SSR; ssr++) {
                hpGraph[inst].gmio_out[ssr].wait();
            }
            //TODO: DEBUG
            print_arr(data_arrays[inst], MAT_ROWS, MAT_COLS);
        }

        for(int inst=0; inst<INSTANCES; inst++) {
            printf("\nPERFORM IFFT (ITER = %d) (INST = %d)\n", iter, inst);
            ifftGraph[inst].run(MAT_ROWS);
            for(int row=0; row<MAT_ROWS; row++) {
                for(int ssr=0; ssr<FFT_NPORTS; ssr++) {
                    ifftGraph[inst].gmio_in[ssr].gm2aie_nb(data_arrays[inst] + MAT_COLS*row + ssr*fft_per_ssr_entry_size, fft_per_ssr_byte_size);
                    ifftGraph[inst].gmio_out[ssr].aie2gm_nb(data_arrays[inst] + MAT_COLS*row + ssr*fft_per_ssr_entry_size, fft_per_ssr_byte_size);
                }
            }
        }

        // Block until AIE finishes
        for(int inst=0; inst<INSTANCES; inst++) {
            for(int ssr=0; ssr<FFT_NPORTS; ssr++) {
                ifftGraph[inst].gmio_out[ssr].wait();
            }
            total_err_cnt+=ifftErrorCheck(data_arrays[inst], gold_data_array);
        }
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
