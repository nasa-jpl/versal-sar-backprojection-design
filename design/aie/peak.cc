// By: Austin Owens
// Date: 6/25/2024
// Desc: Performs a complex conjugate

#include "aie_api/aie.hpp"
#include "aie_api/aie_adf.hpp"
#include "aie_api/utils.hpp"
#include "../common.h"

using namespace adf;


void cols_peak_kern(input_buffer<TT_DATA, extents<MAT_COLS>> &restrict in, 
                    output_buffer<TT_DATA, extents<1>> &restrict out) {
    auto in_iter = aie::cbegin(in);
    auto out_iter = aie::begin(out);

    TT_DATA col_max_val_idx = {0, 0};
    float in_val = 0;
    for(uint16 col = 0; col < MAT_COLS; col++) chess_prepare_for_pipelining {
        //in_val = aie::reduce_max(aie::real(*in_iter));
        in_val = aie::real(*in_iter++);
        if (in_val > col_max_val_idx.real) {
            col_max_val_idx.real = in_val;
            col_max_val_idx.imag = col;
        }
    }

    *out_iter++ = col_max_val_idx;

}

void rows_peak_kern(input_buffer<TT_DATA, extents<MAT_ROWS>> &restrict in, 
                    output_buffer<TT_DATA, extents<1>> &restrict out) {
    auto in_iter = aie::cbegin(in);
    auto out_iter = aie::begin(out);

    TT_DATA max_val_idx = {0, 0};
    TT_DATA col_max_val_idx;
    for(uint16 row = 0; row < MAT_ROWS; row++) chess_prepare_for_pipelining {
        col_max_val_idx = *in_iter++;
        if (col_max_val_idx.real > max_val_idx.real) {
            max_val_idx.real = col_max_val_idx.real;

            // The imag number is actually the index from cols_peak_kern above
            max_val_idx.imag = row*MAT_ROWS + col_max_val_idx.imag;
        }
    }

    *out_iter++ = max_val_idx;

}
