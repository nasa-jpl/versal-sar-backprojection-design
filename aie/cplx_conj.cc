// By: Austin Owens
// Date: 6/25/2024
// Desc: Performs a complex conjugate

#include "aie_api/aie.hpp"
#include "aie_api/aie_adf.hpp"
#include "aie_api/utils.hpp"
#include "../common.h"

using namespace adf;

void cplx_conj_kern(input_buffer<TT_DATA, extents<2048>> &restrict in, 
                    output_buffer<TT_DATA, extents<2048>> &restrict out) {
    auto in_iter = aie::begin_vector<16>(in);
    auto out_iter = aie::begin_vector<16>(out);

    for(unsigned i = 0; i < 2048/16; i++)
		chess_prepare_for_pipelining
    {
        *out_iter++ = aie::conj(*in_iter++);

        //auto data = aie::conj(*in_iter++);
        //printf("i: %d\n", i);
        //aie::print(data, true, "Data:\n");
        //*out_iter++ = data;
    }
}
