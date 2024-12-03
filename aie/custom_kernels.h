// By: Austin Owens
// Date: 6/25/2024
// Desc: Kernel header

#pragma once

#include <adf.h>
#include "aie_api/aie.hpp"
#include "../common.h"

using namespace adf;

void cplx_conj_kern(input_buffer<TT_DATA, extents<2048>>& __restrict in, 
                    output_buffer<TT_DATA, extents<2048>>& __restrict out);

void slowtime_splicer_kern(input_buffer<float, extents<1>>& __restrict x_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict y_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict z_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict ref_range_in,
                           output_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_out);

class ImgReconstruct
{
    public:
        static constexpr float PI = 3.1415926535898;
        static constexpr float TWO_PI = 6.2831853071796;
        static constexpr float INV_TWO_PI = 0.15915494309189;
        static constexpr float C = 299792458.0;

        ImgReconstruct(int id);
        void img_reconstruct_kern(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
                                  input_buffer<TT_DATA, extents<TP_POINT_SIZE/4>>& __restrict xy_px_in,
                                  input_async_buffer<TT_DATA, extents<TP_POINT_SIZE/4>>& __restrict rc_in,
                                  output_async_buffer<TT_DATA, extents<TP_POINT_SIZE/4>>& __restrict img_out,
                                  int32 rtp_pulses_in, int32 &rtp_valid_low_bound_out, int32 &rtp_valid_high_bound_out);

        static void registerKernelClass()
        { 
            REGISTER_FUNCTION(ImgReconstruct::img_reconstruct_kern);
        }

    private:
        uint32 m_id;
        uint32 m_iter;
        cfloat m_prev_low_rc;
        alignas(aie::vector_decl_align) TT_DATA m_img[TP_POINT_SIZE/4];

};

