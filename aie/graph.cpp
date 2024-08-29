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

uint8_t fft_rows_graph_insts = 0;
uint8_t fft_cols_graph_insts = 0;
uint8_t ifft_cols_graph_insts = 0;
uint8_t ifft_rows_graph_insts = 0;
uint8_t cplx_conj_graph_insts = 0;
uint8_t hp_graph_insts = 0;
uint8_t peak_graph_insts = 0;

const int INSTANCES = 4;
FFTRowsGraph fftRowsGraph[INSTANCES];
FFTColsGraph fftColsGraph[INSTANCES];
IFFTColsGraph ifftColsGraph[INSTANCES];
IFFTRowsGraph ifftRowsGraph[INSTANCES];
CplxConjGraph cplxConjGraph[INSTANCES];
HPGraph hpGraph[INSTANCES];
PeakGraph peakGraph[INSTANCES];

#if defined(__AIESIM__) || defined(__X86SIM__)

const int BLOCK_SIZE_ENTRIES = MAT_ROWS * MAT_COLS;
const int BLOCK_SIZE_BYTES = BLOCK_SIZE_ENTRIES * sizeof(TT_DATA);

int fftRowErrorCheck(TT_DATA* fft_row_array, int instance = -1) {
    int error_cnt = 0;

    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
        if(i < MAT_ROWS) {
            if(fft_row_array[i].real != FFT_SAMPLE_DATA || fft_row_array[i].imag != FFT_SAMPLE_DATA)
                error_cnt++;
        } else {
            if(fft_row_array[i].real != 0 || fft_row_array[i].imag != 0)
                error_cnt++;
        }
    }

    if (instance >= 0)
        printf("INSTANCE %d: ", instance);
    printf("FFT ROW ERRORS: %d\n", error_cnt);

    return error_cnt;
}

int fftColErrorCheck(TT_DATA* fft_col_array, int instance = -1) {
    int error_cnt = 0;

    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
        if(fft_col_array[i].real != FFT_SAMPLE_DATA || fft_col_array[i].imag != FFT_SAMPLE_DATA)
            error_cnt++;
    }

    if (instance >= 0)
        printf("INSTANCE %d: ", instance);
    printf("FFT COL ERRORS: %d\n", error_cnt);

    return error_cnt;
}

int cplxConjErrorCheck(TT_DATA* cplx_conj_array, int instance = -1) {
    int error_cnt = 0;

    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
        if(cplx_conj_array[i].real != FFT_SAMPLE_DATA || cplx_conj_array[i].imag != -FFT_SAMPLE_DATA)
            error_cnt++;
    }

    if (instance >= 0)
        printf("INSTANCE %d: ", instance);
    printf("COMPLEX CONJUGATE ERRORS: %d\n", error_cnt);

    return error_cnt;
}

int hpErrorCheck(TT_DATA* hp_array, int instance = -1) {
    int error_cnt = 0;

    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
        if(hp_array[i].real != HP_SAMPLE_DATA || hp_array[i].imag != 0)
            error_cnt++;
    }

    if (instance >= 0)
        printf("INSTANCE %d: ", instance);
    printf("ELEMENT-WISE MATRIX MULTIPLY ERRORS: %d\n", error_cnt);

    return error_cnt;
}

int ifftColErrorCheck(TT_DATA* ifft_col_array, int instance = -1) {
    int error_cnt = 0;

    for(int i = 0; i < BLOCK_SIZE_ENTRIES; i++) {
        // First column should be HP_SAMPLE_DATA*MAT_COLS and everywhere else 0
        if(i%MAT_COLS == 0) {
            if(ifft_col_array[i].real != HP_SAMPLE_DATA*MAT_COLS || ifft_col_array[i].imag != 0)
                error_cnt++;
        } else {
            if(ifft_col_array[i].real != 0 || ifft_col_array[i].imag != 0)
                error_cnt++;
        }
    }

    if (instance >= 0)
        printf("INSTANCE %d: ", instance);
    printf("IFFT COL ERRORS: %d\n", error_cnt);

    return error_cnt;
}

int ifftRowErrorCheck(TT_DATA* ifft_row_array, int instance = -1) {
    int error_cnt = 0;

    if(ifft_row_array[0].real != HP_SAMPLE_DATA*BLOCK_SIZE_ENTRIES || ifft_row_array[0].imag != 0)
        error_cnt++;

    for(int i = 1; i < BLOCK_SIZE_ENTRIES; i++) {
        if(ifft_row_array[i].real != 0 || ifft_row_array[i].imag != 0)
            error_cnt++;
    }

    if (instance >= 0)
        printf("INSTANCE %d: ", instance);
    printf("IFFT ROW ERRORS: %d\n", error_cnt);

    return error_cnt;
}

int peakSearchErrorCheck(TT_DATA* peak_array, int instance = -1) {
    int error_cnt = 0;

    int max_val = peak_array[0].real;
    int index = peak_array[0].imag;

    if(max_val != HP_SAMPLE_DATA*BLOCK_SIZE_ENTRIES || index != 0)
        error_cnt++;

    if (instance >= 0)
        printf("INSTANCE %d: ", instance);
    printf("PEAK SEARCH ERRORS: %d\n", error_cnt);

    return error_cnt;
}


int main(int argc, char ** argv) {

    // Iterations to run. If you want to increment this, make sure there is enough sim data 
    // within the ./aie/aiesim_data/ directory to simulate multiple iteratings. You can do
    // this by matching this number with the ITER_CNT number in the sq_matrix_input_generate.py
    // script and running that script to generate the proper amount of data.
    const int ITER = 1;

    // Total error count
    int total_err_cnt = 0;

    // Initialize all AIE graphs
    for(int inst=0; inst<INSTANCES; inst++) {
        fftRowsGraph[inst].init();
        fftColsGraph[inst].init();
        ifftColsGraph[inst].init();
        ifftRowsGraph[inst].init();
        cplxConjGraph[inst].init();
        hpGraph[inst].init();
        peakGraph[inst].init();
    }
    
    // Create a map_fft_array placeholder so we don't need to keep computing the 
    // map fft every time
    TT_DATA* map_fft_array = (TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES);
    
    // Create and populate the map image with input data that can be validated
    // throughout the pipelines
    TT_DATA* map_array = (TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES);
    map_array[0] = (TT_DATA) {FFT_SAMPLE_DATA, FFT_SAMPLE_DATA};
    for(int j = 1; j < BLOCK_SIZE_ENTRIES; j++) {
        map_array[j] = (TT_DATA) {0, 0};
    }
    
    // Malloc space for template images
    std::vector<TT_DATA*> tmpl_arrays;
    for(int inst=0; inst<INSTANCES; inst++) {
        tmpl_arrays.push_back((TT_DATA*) GMIO::malloc(BLOCK_SIZE_BYTES));
    }
    
    // Number of iteration for the AIE graphs to run
    fftRowsGraph[0].run(MAT_ROWS + (MAT_ROWS*ITER));
    fftColsGraph[0].run(MAT_COLS + (MAT_COLS*ITER));
    ifftColsGraph[0].run(MAT_COLS*ITER);
    ifftRowsGraph[0].run(MAT_ROWS*ITER);
    cplxConjGraph[0].run(MAT_ROWS*ITER);
    hpGraph[0].run((TP_DIM/TP_NUM_FRAMES)*ITER);
    peakGraph[0].run(ITER);
    for(int inst=1; inst<INSTANCES; inst++) {
        fftRowsGraph[inst].run(MAT_ROWS*ITER);
        fftColsGraph[inst].run(MAT_COLS*ITER);
        ifftColsGraph[inst].run(MAT_COLS*ITER);
        ifftRowsGraph[inst].run(MAT_ROWS*ITER);
        cplxConjGraph[inst].run(MAT_ROWS*ITER);
        hpGraph[inst].run((TP_DIM/TP_NUM_FRAMES)*ITER);
        peakGraph[inst].run(ITER);
    }


    std::cout << "\nPERFORM 2D FFT OF MAP IMAGE" << std::endl;
    // Perform row-wise 1D FFT on map image
    fftRowsGraph[0].gmio_in.gm2aie_nb(map_array, BLOCK_SIZE_BYTES);
    fftRowsGraph[0].gmio_out.aie2gm_nb(map_fft_array, BLOCK_SIZE_BYTES);
    fftRowsGraph[0].gmio_out.wait();
    total_err_cnt+=fftRowErrorCheck(map_fft_array);

    // Perform col-wise 1D FFT of map image (PL sends AIE the transposed data)
    fftColsGraph[0].gmio_out.aie2gm_nb(map_fft_array, BLOCK_SIZE_BYTES);
    fftColsGraph[0].gmio_out.wait();
    total_err_cnt+=fftColErrorCheck(map_fft_array);

    // Start image cross correlations between the map and template images for ITER iterations
    for(int n=0; n<ITER; n++) {
        // Create and populate the template images with input data that can be validated
        // throughout the pipeline
        for(int inst=0; inst<INSTANCES; inst++) {
            tmpl_arrays[inst][0] = (TT_DATA) {FFT_SAMPLE_DATA, FFT_SAMPLE_DATA};
            for(int j = 1; j < BLOCK_SIZE_ENTRIES; j++) {
                tmpl_arrays[inst][j] = (TT_DATA) {0, 0};
            }
        }


        std::cout << "\nPERFORM 2D FFT OF TEMPLATE IMAGES" << std::endl;
        // Perform row-wise 1D FFT on all tmplate images at the same time
        for(int inst=0; inst<INSTANCES; inst++) {
            fftRowsGraph[inst].gmio_in.gm2aie_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
            fftRowsGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
        }

        // Block until AIE sends host computed 1D FFT row-wise results
        for(int inst=0; inst<INSTANCES; inst++) {
            fftRowsGraph[inst].gmio_out.wait();
            total_err_cnt+=fftRowErrorCheck(tmpl_arrays[inst], inst);
        }

        for(int inst=0; inst<INSTANCES; inst++) {
        }

        // Perform col-wise 1D FFT on all template images (PL sends AIE the transposed data)
        for(int inst=0; inst<INSTANCES; inst++) {
            fftColsGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
        }

        // Block until AIE finishes
        for(int inst=0; inst<INSTANCES; inst++) {
            fftColsGraph[inst].gmio_out.wait();
            total_err_cnt+=fftColErrorCheck(tmpl_arrays[inst], inst);
        }


        std::cout << "\nPERFORM COMPLEX CONJUGATE ON TEMPLATE IMAGES" << std::endl;
        // Perform complex conjugate on all tmplate images at the same time
        for(int inst=0; inst<INSTANCES; inst++) {
            cplxConjGraph[inst].gmio_in.gm2aie_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
            cplxConjGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
        }

        // Block until AIE finishes
        for(int inst=0; inst<INSTANCES; inst++) {
            cplxConjGraph[inst].gmio_out.wait();
            total_err_cnt+=cplxConjErrorCheck(tmpl_arrays[inst], inst);
        }

        std::cout << "\nPERFORM ELEMENT-WISE MATRIX MULTIPLY BETWEEN MAP AND TEMPLATE IMAGES" << std::endl;
        // Perform element-wise matrix multiply te on all tmplate images at the same time
        int per_ssr_byte_size = BLOCK_SIZE_BYTES / TP_SSR;
        int per_ssr_entry_size = BLOCK_SIZE_ENTRIES / TP_SSR;
        for(int inst=0; inst<INSTANCES; inst++) {
            for(int ssr=0; ssr<TP_SSR; ssr++) {
                hpGraph[inst].gmio_in_A[ssr].gm2aie_nb(map_fft_array, per_ssr_byte_size);
                hpGraph[inst].gmio_in_B[ssr].gm2aie_nb(tmpl_arrays[inst], per_ssr_byte_size);
                hpGraph[inst].gmio_out[ssr].aie2gm_nb(tmpl_arrays[inst]+ssr*per_ssr_entry_size, per_ssr_byte_size);
            }
        }

        // Block until AIE finishes
        for(int inst=0; inst<INSTANCES; inst++) {
            for(int ssr=0; ssr<TP_SSR; ssr++) {
                hpGraph[inst].gmio_out[ssr].wait();
            }
            total_err_cnt+=hpErrorCheck(tmpl_arrays[inst], inst);
        }


        std::cout << "\nPERFORM 2D iFFT OF IMAGES" << std::endl;
        // Perform col-wise 1D iFFT on all images at the same time
        for(int inst=0; inst<INSTANCES; inst++) {
            ifftColsGraph[inst].gmio_in.gm2aie_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
            ifftColsGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
        }

        // Block until AIE finishes
        for(int inst=0; inst<INSTANCES; inst++) {
            ifftColsGraph[inst].gmio_out.wait();
            total_err_cnt+=ifftColErrorCheck(tmpl_arrays[inst], inst);
        }

        // Perform row-wise 1D iFFT on all images (PL sends AIE the transposed data)
        for(int inst=0; inst<INSTANCES; inst++) {
            ifftRowsGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
        }

        // Block until AIE finishes
        for(int inst=0; inst<INSTANCES; inst++) {
            ifftRowsGraph[inst].gmio_out.wait();
            total_err_cnt+=ifftRowErrorCheck(tmpl_arrays[inst], inst);
        }


        std::cout << "\nPERFORM PEAK SEARCH ON ALL IMAGES" << std::endl;
        // Perform peak search on all images at the same time
        for(int inst=0; inst<INSTANCES; inst++) {
            peakGraph[inst].gmio_in.gm2aie_nb(tmpl_arrays[inst], BLOCK_SIZE_BYTES);
            peakGraph[inst].gmio_out.aie2gm_nb(tmpl_arrays[inst], sizeof(TT_DATA));
        }

        // Block until AIE finishes
        for(int inst=0; inst<INSTANCES; inst++) {
            peakGraph[inst].gmio_out.wait();
            total_err_cnt+=peakSearchErrorCheck(tmpl_arrays[inst], inst);
        }

    }

    if(total_err_cnt) {
        printf("\nTEST FAILED!!!\n\n");
    } else {
        printf("\nTEST PASSED!!!\n\n");
    }

    // Free memory
    GMIO::free(map_array);
    GMIO::free(map_fft_array);
    for(int inst=0; inst<INSTANCES; inst++) {
        GMIO::free(tmpl_arrays[inst]);

        // End AIE graphs
        fftRowsGraph[inst].end();
        fftColsGraph[inst].end();
        ifftRowsGraph[inst].end();
        ifftColsGraph[inst].end();
        cplxConjGraph[inst].end();
        hpGraph[inst].end();
        peakGraph[inst].end();
    }

    return total_err_cnt;
}

#endif
