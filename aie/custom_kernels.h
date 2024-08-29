// By: Austin Owens
// Date: 6/25/2024
// Desc: Complex conjugate kernel header

#ifndef CPLXCONJ_H
#define CPLXCONJ_H

#include <adf.h>
#include "../common.h"

using namespace adf;

void cplx_conj_kern(input_buffer<TT_DATA, extents<TP_POINT_SIZE>> & restrict in, 
                    output_buffer<TT_DATA, extents<TP_POINT_SIZE>> & restrict out);

void cols_peak_kern(input_buffer<TT_DATA, extents<MAT_COLS>> & restrict in, 
                    output_buffer<TT_DATA, extents<1>> & restrict out);

void rows_peak_kern(input_buffer<TT_DATA, extents<MAT_ROWS>> & restrict in, 
                    output_buffer<TT_DATA, extents<1>> & restrict out);


#endif
