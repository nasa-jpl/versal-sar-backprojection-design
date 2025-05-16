// By: Austin Owens
// Date: 6/25/2024
// Desc: Kernel header

#pragma once

#include <adf.h>
#include "aie_api/aie.hpp"
#include "../common.h"

using namespace adf;

void cplx_conj_kern(input_buffer<cfloat, extents<2048>>& __restrict in, 
                    output_buffer<cfloat, extents<2048>>& __restrict out);


void px_arbiter_kern(input_stream<float>* __restrict px_xyz_in,
                     output_pktstream *px_xyz_out);

void slowtime_splicer_kern(input_buffer<float, extents<1>>& __restrict x_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict y_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict z_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict ref_range_in,
                           input_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_in,
                           output_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_out,
                           output_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_out);

//void arbiter_kern(input_pktstream *in, output_stream<int> *__restrict out);
//void arbiter_kern(input_pktstream *in, output_pktstream *out);

class ImgReconstruct
{
    public:
        ImgReconstruct(int id);
        void img_reconstruct_kern(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
                                  input_async_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_in,
                                  input_pktstream *px_xyz_in,
                                  //output_pktstream *img_out,
                                  output_async_buffer<cfloat, extents<(PULSES*RC_SAMPLES)/IMG_SOLVERS>>& __restrict img_out,
                                  int rtp_dump_img_in);

        static void registerKernelClass()
        { 
            REGISTER_FUNCTION(ImgReconstruct::img_reconstruct_kern);
        }

    private:
        uint32 m_id;

        // Used for focusing image over several iterations
        alignas(aie::vector_decl_align) cfloat m_img[(PULSES*RC_SAMPLES)/IMG_SOLVERS];
};

//class ImgReconstructA
//{
//    public:
//        static constexpr float PI = 3.1415926535898;
//        static constexpr float TWO_PI = 6.2831853071796;
//        static constexpr float INV_TWO_PI = 0.15915494309189;
//        static constexpr float C = 299792458.0;
//
//        ImgReconstructA(int id);
//        void img_reconstruct_kern_a(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
//                                    input_async_buffer<cfloat, extents<(TP_POINT_SIZE/4)+1>>& __restrict rc_in,
//                                    output_stream<cfloat>* __restrict img_out,
//                                    int32 &rtp_img_elem_cnt_out);
//
//        static void registerKernelClass()
//        { 
//            REGISTER_FUNCTION(ImgReconstructA::img_reconstruct_kern_a);
//        }
//
//    private:
//        uint32 m_id;
//        uint32 m_iter;
//        cfloat m_prev_low_rc;
//        cfloat m_prev_img_val;
//        alignas(aie::vector_decl_align) cfloat m_img[TP_POINT_SIZE/4];
//};
//
//class ImgReconstructB
//{
//    public:
//        static constexpr float PI = 3.1415926535898;
//        static constexpr float TWO_PI = 6.2831853071796;
//        static constexpr float INV_TWO_PI = 0.15915494309189;
//        static constexpr float C = 299792458.0;
//
//        ImgReconstructB(int id);
//        void img_reconstruct_kern_b(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
//                                    input_async_buffer<cfloat, extents<(TP_POINT_SIZE/4)+1>>& __restrict rc_in,
//                                    input_stream<cfloat>* __restrict img_in, 
//                                    output_stream<cfloat>* __restrict img_out,
//                                    int32 &rtp_img_elem_cnt_out);
//
//        static void registerKernelClass()
//        { 
//            REGISTER_FUNCTION(ImgReconstructB::img_reconstruct_kern_b);
//        }
//
//    private:
//        uint32 m_id;
//        uint32 m_iter;
//        cfloat m_prev_low_rc;
//        cfloat m_prev_img_val;
//        alignas(aie::vector_decl_align) cfloat m_img[TP_POINT_SIZE/4];
//};
//
//class ImgReconstructC
//{
//    public:
//        static constexpr float PI = 3.1415926535898;
//        static constexpr float TWO_PI = 6.2831853071796;
//        static constexpr float INV_TWO_PI = 0.15915494309189;
//        static constexpr float C = 299792458.0;
//
//        ImgReconstructC(int id);
//        void img_reconstruct_kern_c(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
//                                    input_async_buffer<cfloat, extents<(TP_POINT_SIZE/4)+1>>& __restrict rc_in,
//                                    input_stream<cfloat>* __restrict img_in, 
//                                    output_stream<cfloat>* __restrict img_out,
//                                    int32 &rtp_img_elem_cnt_out);
//
//        static void registerKernelClass()
//        { 
//            REGISTER_FUNCTION(ImgReconstructC::img_reconstruct_kern_c);
//        }
//
//    private:
//        uint32 m_id;
//        uint32 m_iter;
//        cfloat m_prev_low_rc;
//        cfloat m_prev_img_val;
//        alignas(aie::vector_decl_align) cfloat m_img[TP_POINT_SIZE/4];
//};
