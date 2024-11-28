// By: Austin Owens
// Date: 6/3/2024
// Desc: Using Vitis DSP lib to perform 1D FFT operation

#pragma once

#include "adf.h"
#include "fft_ifft_dit_1ch_graph.hpp"
#include "hadamard_graph.hpp"
#include "custom_kernels.h"

using namespace adf;
namespace dsplib = xf::dsp::aie;

extern uint8_t fft_graph_insts;
extern uint8_t ifft_graph_insts;
extern uint8_t hp_graph_insts;
extern uint8_t cplx_conj_graph_insts;
extern uint8_t bp_graph_insts;

//template<int FFT_X, int FFT_Y>
//class FFTGraph: public graph
//{
//    public:
//        input_gmio gmio_in[FFT_NPORTS];
//        output_gmio gmio_out[FFT_NPORTS];
//
//      
//        FFTGraph() {
//            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
//                                                         TT_TWIDDLE,
//                                                         TP_POINT_SIZE,
//                                                         TP_FFT,
//                                                         TP_FFT_SHIFT,
//                                                         TP_FFT_CASC_LEN,
//                                                         TP_DYN_PT_SIZE,
//                                                         TP_FFT_WINDOW_VSIZE,
//                                                         TP_FFT_API,
//                                                         TP_PARALLEL_POWER,
//                                                         TP_USE_WIDGETS> fftGraph;
//
//            // Runtime ratio gives the tools the flexibility to put multiple AIE kernels on a single tile
//            runtime<ratio>(*fftGraph.FFTsubframe0.FFTsubframe0.getKernels()) = 1.0;
//            runtime<ratio>(*fftGraph.FFTsubframe0.FFTsubframe1.getKernels()) = 1.0;
//            runtime<ratio>(*fftGraph.FFTsubframe1.FFTsubframe0.getKernels()) = 1.0;
//            runtime<ratio>(*fftGraph.FFTsubframe1.FFTsubframe1.getKernels()) = 1.0;
//
//            for (int i = 0; i < FFT_NPORTS; i++) {
//                gmio_in[i] = input_gmio::create("gmio_in_fft_" + std::to_string(fft_graph_insts) + "_" + std::to_string(i), 256, 1000);
//                connect< window<FFT_WINDOW_BUFF_SIZE> >(gmio_in[i].out[0],  fftGraph.in[i]);
//
//                gmio_out[i] = output_gmio::create("gmio_out_fft_" + std::to_string(fft_graph_insts) + "_" + std::to_string(i), 256, 1000);
//                connect< window<FFT_WINDOW_BUFF_SIZE> >(fftGraph.out[i], gmio_out[i].in[0]);
//            }
//
//            if (FFT_X > 0 && FFT_Y > 0) {
//                if (TP_FFT_CASC_LEN == 4) {
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[0]) = tile(FFT_X + 0,  FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[1]) = tile(FFT_X + 1,  FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[2]) = tile(FFT_X + 2,  FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[3]) = tile(FFT_X + 3,  FFT_Y - 0);
//
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[3]) = tile(FFT_X + 6,  FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[2]) = tile(FFT_X + 7,  FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[1]) = tile(FFT_X + 8,  FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[0]) = tile(FFT_X + 9,  FFT_Y - 0);
//
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[0]) = tile(FFT_X + 10, FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[1]) = tile(FFT_X + 11, FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[2]) = tile(FFT_X + 12, FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[3]) = tile(FFT_X + 13, FFT_Y - 0);
//
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[3]) = tile(FFT_X + 16, FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[2]) = tile(FFT_X + 17, FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[1]) = tile(FFT_X + 18, FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[0]) = tile(FFT_X + 19, FFT_Y - 0);
//
//                    location<kernel>(fftGraph.FFTsubframe0.m_r2Comb[0]) = tile(FFT_X + 4,  FFT_Y - 1);
//                    location<kernel>(fftGraph.FFTsubframe0.m_r2Comb[1]) = tile(FFT_X + 5,  FFT_Y - 1);
//
//                    location<kernel>(fftGraph.FFTsubframe1.m_r2Comb[0]) = tile(FFT_X + 14, FFT_Y - 1);
//                    location<kernel>(fftGraph.FFTsubframe1.m_r2Comb[1]) = tile(FFT_X + 15, FFT_Y - 1);
//
//                    location<kernel>(fftGraph.m_r2Comb[0]) = tile(FFT_X + 8,  FFT_Y - 2);
//                    location<kernel>(fftGraph.m_r2Comb[1]) = tile(FFT_X + 9,  FFT_Y - 2);
//                    location<kernel>(fftGraph.m_r2Comb[2]) = tile(FFT_X + 10, FFT_Y - 2);
//                    location<kernel>(fftGraph.m_r2Comb[3]) = tile(FFT_X + 11, FFT_Y - 2);
//                } else if (TP_FFT_CASC_LEN == 1) {
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[0]) = tile(FFT_X + 0, FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[0]) = tile(FFT_X + 3, FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[0]) = tile(FFT_X + 6, FFT_Y - 0);
//                    location<kernel>(fftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[0]) = tile(FFT_X + 9, FFT_Y - 0);
//
//                    location<kernel>(fftGraph.FFTsubframe0.m_r2Comb[0]) = tile(FFT_X + 1,  FFT_Y - 1);
//                    location<kernel>(fftGraph.FFTsubframe0.m_r2Comb[1]) = tile(FFT_X + 2,  FFT_Y - 1);
//
//                    location<kernel>(fftGraph.FFTsubframe1.m_r2Comb[0]) = tile(FFT_X + 7, FFT_Y - 1);
//                    location<kernel>(fftGraph.FFTsubframe1.m_r2Comb[1]) = tile(FFT_X + 8, FFT_Y - 1);
//
//                    location<kernel>(fftGraph.m_r2Comb[0]) = tile(FFT_X + 3, FFT_Y - 2);
//                    location<kernel>(fftGraph.m_r2Comb[1]) = tile(FFT_X + 4, FFT_Y - 2);
//                    location<kernel>(fftGraph.m_r2Comb[2]) = tile(FFT_X + 5, FFT_Y - 2);
//                    location<kernel>(fftGraph.m_r2Comb[3]) = tile(FFT_X + 6, FFT_Y - 2);
//                }
//            }
//
//            ++fft_graph_insts;
//        }
//};
//
//template<int IFFT_X, int IFFT_Y>
//class IFFTGraph: public graph
//{
//    public:
//        input_gmio gmio_in[FFT_NPORTS];
//        output_gmio gmio_out[FFT_NPORTS];
//        
//        IFFTGraph() {
//            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
//                                                         TT_TWIDDLE,
//                                                         TP_POINT_SIZE,
//                                                         TP_IFFT,
//                                                         TP_FFT_SHIFT,
//                                                         TP_FFT_CASC_LEN,
//                                                         TP_DYN_PT_SIZE,
//                                                         TP_FFT_WINDOW_VSIZE,
//                                                         TP_FFT_API,
//                                                         TP_PARALLEL_POWER,
//                                                         TP_USE_WIDGETS> ifftGraph;
//            
//            // Runtime ratio gives the tools the flexibility to put multiple AIE kernels on a single tile
//            runtime<ratio>(*ifftGraph.FFTsubframe0.FFTsubframe0.getKernels()) = 1.0;
//            runtime<ratio>(*ifftGraph.FFTsubframe0.FFTsubframe1.getKernels()) = 1.0;
//            runtime<ratio>(*ifftGraph.FFTsubframe1.FFTsubframe0.getKernels()) = 1.0;
//            runtime<ratio>(*ifftGraph.FFTsubframe1.FFTsubframe1.getKernels()) = 1.0;
//    
//            for (int i = 0; i < FFT_NPORTS; i++) {
//                gmio_in[i] = input_gmio::create("gmio_in_ifft_" + std::to_string(ifft_graph_insts) + "_" + std::to_string(i), 256, 1000);
//                connect< window<FFT_WINDOW_BUFF_SIZE> >(gmio_in[i].out[0], ifftGraph.in[i]);
//
//                gmio_out[i] = output_gmio::create("gmio_out_ifft_" + std::to_string(ifft_graph_insts) + "_" + std::to_string(i), 256, 1000);
//                connect< window<FFT_WINDOW_BUFF_SIZE> >(ifftGraph.out[i], gmio_out[i].in[0]);
//            }
//            
//            if (IFFT_X > 0 && IFFT_Y > 0) {
//                if (TP_FFT_CASC_LEN == 4) {
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[0]) = tile(IFFT_X + 0,  IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[1]) = tile(IFFT_X + 1,  IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[2]) = tile(IFFT_X + 2,  IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[3]) = tile(IFFT_X + 3,  IFFT_Y - 0);
//
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[3]) = tile(IFFT_X + 6,  IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[2]) = tile(IFFT_X + 7,  IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[1]) = tile(IFFT_X + 8,  IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[0]) = tile(IFFT_X + 9,  IFFT_Y - 0);
//
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[0]) = tile(IFFT_X + 10, IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[1]) = tile(IFFT_X + 11, IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[2]) = tile(IFFT_X + 12, IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[3]) = tile(IFFT_X + 13, IFFT_Y - 0);
//
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[3]) = tile(IFFT_X + 16, IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[2]) = tile(IFFT_X + 17, IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[1]) = tile(IFFT_X + 18, IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[0]) = tile(IFFT_X + 19, IFFT_Y - 0);
//
//                    location<kernel>(ifftGraph.FFTsubframe0.m_r2Comb[0]) = tile(IFFT_X + 4,  IFFT_Y - 1);
//                    location<kernel>(ifftGraph.FFTsubframe0.m_r2Comb[1]) = tile(IFFT_X + 5,  IFFT_Y - 1);
//
//                    location<kernel>(ifftGraph.FFTsubframe1.m_r2Comb[0]) = tile(IFFT_X + 14, IFFT_Y - 1);
//                    location<kernel>(ifftGraph.FFTsubframe1.m_r2Comb[1]) = tile(IFFT_X + 15, IFFT_Y - 1);
//
//                    location<kernel>(ifftGraph.m_r2Comb[0]) = tile(IFFT_X + 8,  IFFT_Y - 2);
//                    location<kernel>(ifftGraph.m_r2Comb[1]) = tile(IFFT_X + 9,  IFFT_Y - 2);
//                    location<kernel>(ifftGraph.m_r2Comb[2]) = tile(IFFT_X + 10, IFFT_Y - 2);
//                    location<kernel>(ifftGraph.m_r2Comb[3]) = tile(IFFT_X + 11, IFFT_Y - 2);
//                } 
//                else if (TP_FFT_CASC_LEN == 1) {
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe0.FFTwinproc.m_fftKernels[0]) = tile(IFFT_X + 0, IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe1.FFTwinproc.m_fftKernels[0]) = tile(IFFT_X + 3, IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe0.FFTwinproc.m_fftKernels[0]) = tile(IFFT_X + 6, IFFT_Y - 0);
//                    location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe1.FFTwinproc.m_fftKernels[0]) = tile(IFFT_X + 9, IFFT_Y - 0);
//
//                    location<kernel>(ifftGraph.FFTsubframe0.m_r2Comb[0]) = tile(IFFT_X + 1,  IFFT_Y - 1);
//                    location<kernel>(ifftGraph.FFTsubframe0.m_r2Comb[1]) = tile(IFFT_X + 2,  IFFT_Y - 1);
//
//                    location<kernel>(ifftGraph.FFTsubframe1.m_r2Comb[0]) = tile(IFFT_X + 7, IFFT_Y - 1);
//                    location<kernel>(ifftGraph.FFTsubframe1.m_r2Comb[1]) = tile(IFFT_X + 8, IFFT_Y - 1);
//
//                    location<kernel>(ifftGraph.m_r2Comb[0]) = tile(IFFT_X + 3, IFFT_Y - 2);
//                    location<kernel>(ifftGraph.m_r2Comb[1]) = tile(IFFT_X + 4, IFFT_Y - 2);
//                    location<kernel>(ifftGraph.m_r2Comb[2]) = tile(IFFT_X + 5, IFFT_Y - 2);
//                    location<kernel>(ifftGraph.m_r2Comb[3]) = tile(IFFT_X + 6, IFFT_Y - 2);
//                }
//            }
//    
//            ++ifft_graph_insts;
//        }
//};
//
//class HPGraph: public graph
//{
//    public:
//        input_gmio gmio_in_A[TP_SSR];
//        input_gmio gmio_in_B[TP_SSR];
//        output_gmio gmio_out[TP_SSR];
//      
//        HPGraph() {
//            dsplib::hadamard::hadamard_graph<TT_DATA_A,
//                                             TT_DATA_B,
//                                             TP_DIM,
//                                             TP_NUM_FRAMES,
//                                             TP_HP_SHIFT,
//                                             TP_HP_API,
//                                             TP_SSR,
//                                             TP_RND,
//                                             TP_SAT> hpGraph;
//         
//            // Runtime ratio gives the tools the flexibility to put multiple AIE kernels on a single tile
//            runtime<ratio>(*hpGraph.getKernels()) = 1.0;
//
//            for (int i = 0; i < TP_SSR; i++) {
//                gmio_in_A[i] = input_gmio::create("gmio_in_vecA_" + std::to_string(hp_graph_insts) + "_" + std::to_string(i), 256, 1000);
//                gmio_in_B[i] = input_gmio::create("gmio_in_vecB_" + std::to_string(hp_graph_insts) + "_" + std::to_string(i), 256, 1000);
//
//                connect(gmio_in_A[i].out[0], hpGraph.inA[i]);
//                connect(gmio_in_B[i].out[0], hpGraph.inB[i]);
//
//                gmio_out[i] = output_gmio::create("gmio_out_vec_" + std::to_string(hp_graph_insts) + "_" + std::to_string(i), 256, 1000);
//                connect(hpGraph.out[i], gmio_out[i].in[0]);
//            }
//
//            ++hp_graph_insts;
//        }
//};
//
//class CplxConjGraph: public graph
//{
//    private:
//        kernel k_m;
//
//    public:
//        output_gmio gmio_out;
//        input_gmio gmio_in;
//
//        CplxConjGraph() {
//            k_m = kernel::create(cplx_conj_kern);
//            gmio_out = output_gmio::create("gmio_out_cplx_conj_" + std::to_string(cplx_conj_graph_insts), 256, 1000);
//            gmio_in = input_gmio::create("gmio_in_cplx_conj_" + std::to_string(cplx_conj_graph_insts), 256, 1000);
//
//            connect(gmio_in.out[0], k_m.in[0]);
//            connect(k_m.out[0], gmio_out.in[0]);
//
//            source(k_m) = "cplx_conj.cc";
//
//            // Runtime ratio gives the tools the flexibility to put multiple AIE kernels on a single tile
//            runtime<ratio>(k_m) = 1.0;
//
//            ++cplx_conj_graph_insts;
//        }
//};
//
////class PhaseCorrGraph: public graph
////{
////    private:
////        kernel k_m;
////
////    public:
////        input_gmio gmio_in_x_ant_pos;
////        input_gmio gmio_in_y_ant_pos;
////        input_gmio gmio_in_z_ant_pos;
////        input_gmio gmio_in_ref_range;
////        input_gmio gmio_in_ph_data;
////        output_gmio gmio_out_img;
////
////        PhaseCorrGraph() {
////            k_m = kernel::create(phase_corr_kern);
////            gmio_in_x_ant_pos = input_gmio::create("gmio_in_x_ant_pos_" + std::to_string(phase_corr_graph_insts), 256, 1000);
////            gmio_in_y_ant_pos = input_gmio::create("gmio_in_y_ant_pos_" + std::to_string(phase_corr_graph_insts), 256, 1000);
////            gmio_in_z_ant_pos = input_gmio::create("gmio_in_z_ant_pos_" + std::to_string(phase_corr_graph_insts), 256, 1000);
////            gmio_in_ref_range = input_gmio::create("gmio_in_ref_range_" + std::to_string(phase_corr_graph_insts), 256, 1000);
////            gmio_in_ph_data = input_gmio::create("gmio_in_ph_data_" + std::to_string(phase_corr_graph_insts), 256, 1000);
////            gmio_out_img = output_gmio::create("gmio_out_img_" + std::to_string(phase_corr_graph_insts), 256, 1000);
////
////            connect<>(gmio_in_x_ant_pos.out[0], k_m.in[0]);
////            connect<>(gmio_in_y_ant_pos.out[0], k_m.in[1]);
////            connect<>(gmio_in_z_ant_pos.out[0], k_m.in[2]);
////            connect<>(gmio_in_ref_range.out[0], k_m.in[3]);
////            connect<>(gmio_in_ph_data.out[0], k_m.in[4]);
////            connect<>(k_m.out[0], gmio_out_img.in[0]);
////
////            source(k_m) = "phase_corr.cc";
////
////            // Runtime ratio gives the tools the flexibility to put multiple AIE kernels on a single tile
////            runtime<ratio>(k_m) = 1.0;
////
////            ++phase_corr_graph_insts;
////        }
////};

class BackProjectionGraph: public graph
{
    private:

        //***** KERNEL OBJECTS *****//

        // Slow time splicer kernel module
        kernel sts_km;

        // Image reconstruction kernel module
        kernel img_rec_km[IMG_SOLVERS];

    public:
        //***** GMIO PORT OBJECTS *****//

        // Slow time splicer GMIO ports
        input_gmio gmio_in_x_ant_pos;
        input_gmio gmio_in_y_ant_pos;
        input_gmio gmio_in_z_ant_pos;
        input_gmio gmio_in_ref_range;

        // Image reconstruction GMIO ports
        input_gmio gmio_in_xy_px[IMG_SOLVERS];
        input_gmio gmio_in_rc[IMG_SOLVERS];
        output_gmio gmio_out_img[IMG_SOLVERS];
            

        //***** RTP PORT OBJECTS *****//
        inout_port rtp_valid_low_bound[4];
        inout_port rtp_valid_high_bound[4];

        BackProjectionGraph() {

            //***** KERNELS *****//

            // Slow time splicer kernel
            sts_km = kernel::create(slowtime_splicer_kern);

            // Instantiating ImgReconstruct allows for keeping track of instances across 
            // instantiations. This gives each backprojection kernel instance an ID it 
            // can reference to do something unique based on its ID.
            for (int i=0; i<IMG_SOLVERS; i++) {
                img_rec_km[i] = kernel::create_object<ImgReconstruct>(i);
            }


            //***** GMIO PORTS *****//
            
            // Slow time splicer GMIO ports
            gmio_in_x_ant_pos = input_gmio::create("gmio_in_x_ant_pos_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_y_ant_pos = input_gmio::create("gmio_in_y_ant_pos_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_z_ant_pos = input_gmio::create("gmio_in_z_ant_pos_" + std::to_string(bp_graph_insts), 256, 1000);
            gmio_in_ref_range = input_gmio::create("gmio_in_ref_range_" + std::to_string(bp_graph_insts), 256, 1000);

            // Image Reconstruct GMIO ports
            for (int i=0; i<IMG_SOLVERS; i++) {
                gmio_in_xy_px[i] = input_gmio::create("gmio_in_xy_px_" + std::to_string(bp_graph_insts) + "_" + std::to_string(i), 256, 1000);
                gmio_in_rc[i] = input_gmio::create("gmio_in_rc_" + std::to_string(bp_graph_insts) + "_" + std::to_string(i), 256, 1000);
                gmio_out_img[i] = output_gmio::create("gmio_out_img_" + std::to_string(bp_graph_insts) + "_" + std::to_string(i), 256, 1000);
            }


            //***** GMIO CONNECTIONS *****//

            // Slow time splicer GMIO connections
            connect(gmio_in_x_ant_pos.out[0], sts_km.in[0]);
            connect(gmio_in_y_ant_pos.out[0], sts_km.in[1]);
            connect(gmio_in_z_ant_pos.out[0], sts_km.in[2]);
            connect(gmio_in_ref_range.out[0], sts_km.in[3]);


            // Image reconstruct GMIO connections
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect(gmio_in_xy_px[i].out[0], img_rec_km[i].in[1]);
                connect(gmio_in_rc[i].out[0], img_rec_km[i].in[2]);
                connect(img_rec_km[i].out[0], gmio_out_img[i].in[0]);
            }

            //for (int i=0; i<IMG_SOLVERS; i++) {
            //    int fftshift_id = i-(IMG_SOLVERS/2);
            //    if (fftshift_id < 0)
            //        fftshift_id+=IMG_SOLVERS;
            //    connect(gmio_in_rc[i].out[0], img_rec_km[fftshift_id].in[2]);
            //}


            //***** AIE TO AIE CONNECTIONS *****//

            // Slow time splicer to image reconstruction
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect(sts_km.out[0], img_rec_km[i].in[0]);
            }

            //***** RTP CONNECTIONS *****//
            
            // image reconstruction to valid_bounds RTP param
            for (int i=0; i<IMG_SOLVERS; i++) {
                connect<parameter>(sync(img_rec_km[i].inout[0]), rtp_valid_low_bound[i]);
                connect<parameter>(sync(img_rec_km[i].inout[1]), rtp_valid_high_bound[i]);
            }
    

            //***** SOURCE FILES *****//

            source(sts_km) = "backprojection.cc";
            for (int i=0; i<IMG_SOLVERS; i++)
                source(img_rec_km[i]) = "backprojection.cc";


            //***** RUNTIME RATIOS *****//

            runtime<ratio>(sts_km) = 1.0;
            for (int i=0; i<IMG_SOLVERS; i++)
                runtime<ratio>(img_rec_km[i]) = 1.0;


            ++bp_graph_insts;

            //for (int i=0; i<BP_SOLVERS; i++) {
            //    connect(gmio_in_rc[i].out[0], img_rec_km[(BP_SOLVERS-1)-i].in[1]);
            //}
            //for (int i=0; i<BP_SOLVERS; i++) {
            //    // ImgReconstruct GMIO connections accounting for fftshifting
            //    int bp_fftshift_id = i-(BP_SOLVERS/2);
            //    if (bp_fftshift_id < 0)
            //        bp_fftshift_id+=BP_SOLVERS;
            //    connect(gmio_in_rc[i].out[0], img_rec_km[bp_fftshift_id].in[1]);
            //    //connect(gmio_in_xy_px[i].out[0], img_rec_km[bp_fftshift_id].in[2]);
            //    //connect(img_rec_km[bp_fftshift_id].out[0], gmio_out_img[bp_fftshift_id].in[0]);
            //}
        
        }
};
