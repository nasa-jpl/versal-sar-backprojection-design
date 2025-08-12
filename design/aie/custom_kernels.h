// By: Austin Owens
// Date: 6/25/2024
// Desc: Kernel header

#pragma once

#include <adf.h>
#include "aie_api/aie.hpp"
#include "../common.h"

using namespace adf;

void px_demux_kern(input_stream<float>* __restrict px_xyz_in,
                   output_pktstream *px_xyz_out);

void data_broadcast_kern(input_stream<float>* __restrict slowtime_in,
                         input_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_in,
                         output_stream<float>* __restrict slowtime_out,
                         output_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_out);

class ImgReconstruct
{
    public:
        ImgReconstruct(int id);
        void img_reconstruct_kern(input_buffer<float, extents<BC_ELEMENTS>>& __restrict slowtime_in,
                                  input_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_in,
                                  input_pktstream *px_xyz_in,
                                  output_pktstream *img_out,
                                  int rtp_dump_img_in);

        static void registerKernelClass()
        { 
            REGISTER_FUNCTION(ImgReconstruct::img_reconstruct_kern);
        }

    private:
        uint32 m_id;

        // Used for focusing image over several iterations
        alignas(aie::vector_decl_align) cfloat m_img[(AZ_SAMPLES*RC_SAMPLES)/IMG_SOLVERS];
};
