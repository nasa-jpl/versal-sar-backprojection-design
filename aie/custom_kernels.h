// By: Austin Owens
// Date: 6/25/2024
// Desc: Complex conjugate kernel header

#ifndef CPLXCONJ_H
#define CPLXCONJ_H

#include <adf.h>
#include "../common.h"

using namespace adf;

void cplx_conj_kern(input_buffer<TT_DATA, extents<2048>> & restrict in, 
                    output_buffer<TT_DATA, extents<2048>> & restrict out);

#endif
