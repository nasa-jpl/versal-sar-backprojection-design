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
//extern uint8_t cplx_conj_graph_insts;

template<int FFT_X, int FFT_Y>
class FFTGraph: public graph
{
    public:
        //input_gmio gmio_in;
        //output_gmio gmio_out;
        input_gmio gmio_in[FFT_NPORTS];
        output_gmio gmio_out[FFT_NPORTS];

      
        FFTGraph() {
            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
                                                         TT_TWIDDLE,
                                                         TP_POINT_SIZE,
                                                         TP_FFT,
                                                         TP_FFT_SHIFT,
                                                         TP_FFT_CASC_LEN,
                                                         TP_DYN_PT_SIZE,
                                                         TP_FFT_WINDOW_VSIZE,
                                                         TP_FFT_API,
                                                         TP_PARALLEL_POWER,
                                                         TP_USE_WIDGETS> fftGraph;
            //dsplib::fft::mixed_radix_fft::mixed_radix_fft_graph<TT_DATA,
            //                                                    TT_TWIDDLE,
            //                                                    TP_POINT_SIZE,
            //                                                    TP_FFT,
            //                                                    TP_FFT_SHIFT,
            //                                                    TP_FFT_RND,
            //                                                    TP_FFT_SAT,
            //                                                    TP_FFT_WINDOW_VSIZE,
            //                                                    TP_FFT_CASC_LEN> fftGraph;

            //runtime<ratio>(*fftGraph.getKernels()) = 0.8;

            //// These are port names in to and out of the AI Engine
            //gmio_in = input_gmio::create("gmio_in_fft_" + std::to_string(fft_graph_insts), 256, 1000);
            //gmio_out = output_gmio::create("gmio_out_fft_" + std::to_string(fft_graph_insts), 256, 1000);

            //connect< window<FFT_WINDOW_BUFF_SIZE> > (gmio_in.out[0], fftGraph.in[0]);
            //connect< window<FFT_WINDOW_BUFF_SIZE> > (fftGraph.out[0], gmio_out.in[0]);

            for (int i = 0; i < FFT_NPORTS; i++) {
                gmio_in[i] = input_gmio::create("gmio_in_fft_" + std::to_string(fft_graph_insts) + "_" + std::to_string(i), 256, 1000);
                connect< window<FFT_WINDOW_BUFF_SIZE> >(gmio_in[i].out[0],  fftGraph.in[i]);

                gmio_out[i] = output_gmio::create("gmio_out_fft_" + std::to_string(fft_graph_insts) + "_" + std::to_string(i), 256, 1000);
                connect< window<FFT_WINDOW_BUFF_SIZE> >(fftGraph.out[i], gmio_out[i].in[0]);
            }

            ////location<kernel>(fftGraph.FFTsubframe0.FFTwinproc.m_inWidgetKernel)  = tile(FFT_X + 8, FFT_Y + 0 );
            //location<kernel>(fftGraph.FFTsubframe0.FFTwinproc.m_fftKernels[0])   = tile(FFT_X + 7, FFT_Y + 0 );
            //location<kernel>(fftGraph.FFTsubframe0.FFTwinproc.m_fftKernels[1])   = tile(FFT_X + 6, FFT_Y + 0 );
            //location<kernel>(fftGraph.FFTsubframe0.FFTwinproc.m_fftKernels[2])   = tile(FFT_X + 5, FFT_Y + 0 );
            ////location<kernel>(fftGraph.FFTsubframe0.FFTwinproc.m_outWidgetKernel) = tile(FFT_X + 4, FFT_Y + 0 );
            //location<kernel>(fftGraph.m_combInKernel[0])                         = tile(FFT_X + 3, FFT_Y + 0 );
            //location<kernel>(fftGraph.m_r2Comb[0])                               = tile(FFT_X + 2, FFT_Y + 0 );
            //location<kernel>(fftGraph.m_combOutKernel[0])                        = tile(FFT_X + 1, FFT_Y + 0 );
            ////location<kernel>(fftGraph.FFTsubframe1.FFTwinproc.m_inWidgetKernel)  = tile(FFT_X + 0, FFT_Y + 1 );
            //location<kernel>(fftGraph.FFTsubframe1.FFTwinproc.m_fftKernels[0])   = tile(FFT_X + 1, FFT_Y + 1 );
            //location<kernel>(fftGraph.FFTsubframe1.FFTwinproc.m_fftKernels[1])   = tile(FFT_X + 2, FFT_Y + 1 );
            //location<kernel>(fftGraph.FFTsubframe1.FFTwinproc.m_fftKernels[2])   = tile(FFT_X + 3, FFT_Y + 1 );
            ////location<kernel>(fftGraph.FFTsubframe1.FFTwinproc.m_outWidgetKernel) = tile(FFT_X + 4, FFT_Y + 1 );
            //location<kernel>(fftGraph.m_combInKernel[1])                         = tile(FFT_X + 5, FFT_Y + 1 );
            //location<kernel>(fftGraph.m_r2Comb[1])                               = tile(FFT_X + 6, FFT_Y + 1 );
            //location<kernel>(fftGraph.m_combOutKernel[1])                        = tile(FFT_X + 7, FFT_Y + 1 );

            ++fft_graph_insts;
        }
};

template<int IFFT_X, int IFFT_Y>
class IFFTGraph: public graph
{
    public:
        //input_gmio gmio_in;
        //output_gmio gmio_out;
        input_gmio gmio_in[FFT_NPORTS];
        output_gmio gmio_out[FFT_NPORTS];
        
        IFFTGraph() {
            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
                                                         TT_TWIDDLE,
                                                         TP_POINT_SIZE,
                                                         TP_IFFT,
                                                         TP_FFT_SHIFT,
                                                         TP_FFT_CASC_LEN,
                                                         TP_DYN_PT_SIZE,
                                                         TP_FFT_WINDOW_VSIZE,
                                                         TP_FFT_API,
                                                         TP_PARALLEL_POWER,
                                                         TP_USE_WIDGETS> ifftGraph;
            //dsplib::fft::mixed_radix_fft::mixed_radix_fft_graph<TT_DATA,
            //                                                    TT_TWIDDLE,
            //                                                    TP_POINT_SIZE,
            //                                                    TP_IFFT,
            //                                                    TP_FFT_SHIFT,
            //                                                    TP_FFT_RND,
            //                                                    TP_FFT_SAT,
            //                                                    TP_FFT_WINDOW_VSIZE,
            //                                                    TP_FFT_CASC_LEN> ifftGraph;

            //runtime<ratio>(*ifftGraph.getKernels()) = 0.8;
    
            //// These are port names in to and out of the AI Engine
            //gmio_in = input_gmio::create("gmio_in_ifft_" + std::to_string(ifft_graph_insts), 256, 1000);
            //gmio_out = output_gmio::create("gmio_out_ifft_" + std::to_string(ifft_graph_insts), 256, 1000);
    
            //connect< window<FFT_WINDOW_BUFF_SIZE> > (gmio_in.out[0], ifftGraph.in[0]);
            //connect< window<FFT_WINDOW_BUFF_SIZE> > (ifftGraph.out[0], gmio_out.in[0]);

            for (int i = 0; i < FFT_NPORTS; i++) {
                gmio_in[i] = input_gmio::create("gmio_in_ifft_" + std::to_string(ifft_graph_insts) + "_" + std::to_string(i), 256, 1000);
                connect< window<FFT_WINDOW_BUFF_SIZE> >(gmio_in[i].out[0], ifftGraph.in[i]);

                gmio_out[i] = output_gmio::create("gmio_out_ifft_" + std::to_string(ifft_graph_insts) + "_" + std::to_string(i), 256, 1000);
                connect< window<FFT_WINDOW_BUFF_SIZE> >(ifftGraph.out[i], gmio_out[i].in[0]);
            }
            
            // Location placement of AIE kernels
            //for (int lane=0; lane<FFT_NPORTS; lane++) {
            //    location<kernel>(ifftGraph.m_r2Comb[lane]) = tile(IFFT_X + lane * 2, IFFT_Y + TP_FFT_CASC_LEN + 4);
            //}

            //for (int lane=0; lane<FFT_NPORTS/2; lane++) {
            //        location<kernel>(ifftGraph.FFTsubframe0.m_r2Comb[lane]) =
            //            tile(IFFT_X + lane * 2, IFFT_Y + TP_FFT_CASC_LEN + 3);
            //        location<kernel>(ifftGraph.FFTsubframe1.m_r2Comb[lane]) =
            //            tile(IFFT_X + lane * 2 + 16, IFFT_Y + TP_FFT_CASC_LEN + 3);
            //}

            //for (int lane=0; lane<FFT_NPORTS/4; lane++) {
            //        location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe0.m_r2Comb[lane]) =
            //            tile(IFFT_X + lane * 2, IFFT_Y + TP_FFT_CASC_LEN + 2);
            //        location<kernel>(ifftGraph.FFTsubframe0.FFTsubframe1.m_r2Comb[lane]) =
            //            tile(IFFT_X + lane * 2 + 8, IFFT_Y + TP_FFT_CASC_LEN + 2);
            //        location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe0.m_r2Comb[lane]) =
            //            tile(IFFT_X + lane * 2 + 16, IFFT_Y + TP_FFT_CASC_LEN + 2);
            //        location<kernel>(ifftGraph.FFTsubframe1.FFTsubframe1.m_r2Comb[lane]) =
            //            tile(IFFT_X + lane * 2 + 24, IFFT_Y + TP_FFT_CASC_LEN + 2);
            //}

            ////location<kernel>(ifftGraph.FFTsubframe0.FFTwinproc.m_inWidgetKernel)  = tile(IFFT_X + 8, IFFT_Y + 0 );
            //location<kernel>(ifftGraph.FFTsubframe0.FFTwinproc.m_fftKernels[0])   = tile(IFFT_X + 7, IFFT_Y + 0 );
            //location<kernel>(ifftGraph.FFTsubframe0.FFTwinproc.m_fftKernels[1])   = tile(IFFT_X + 6, IFFT_Y + 0 );
            //location<kernel>(ifftGraph.FFTsubframe0.FFTwinproc.m_fftKernels[2])   = tile(IFFT_X + 5, IFFT_Y + 0 );
            ////location<kernel>(ifftGraph.FFTsubframe0.FFTwinproc.m_outWidgetKernel) = tile(IFFT_X + 4, IFFT_Y + 0 );
            //location<kernel>(ifftGraph.m_combInKernel[0])                         = tile(IFFT_X + 3, IFFT_Y + 0 );
            //location<kernel>(ifftGraph.m_r2Comb[0])                               = tile(IFFT_X + 2, IFFT_Y + 0 );
            //location<kernel>(ifftGraph.m_combOutKernel[0])                        = tile(IFFT_X + 1, IFFT_Y + 0 );
            ////location<kernel>(ifftGraph.FFTsubframe1.FFTwinproc.m_inWidgetKernel)  = tile(IFFT_X + 0, IFFT_Y + 1 );
            //location<kernel>(ifftGraph.FFTsubframe1.FFTwinproc.m_fftKernels[0])   = tile(IFFT_X + 1, IFFT_Y + 1 );
            //location<kernel>(ifftGraph.FFTsubframe1.FFTwinproc.m_fftKernels[1])   = tile(IFFT_X + 2, IFFT_Y + 1 );
            //location<kernel>(ifftGraph.FFTsubframe1.FFTwinproc.m_fftKernels[2])   = tile(IFFT_X + 3, IFFT_Y + 1 );
            ////location<kernel>(ifftGraph.FFTsubframe1.FFTwinproc.m_outWidgetKernel) = tile(IFFT_X + 4, IFFT_Y + 1 );
            //location<kernel>(ifftGraph.m_combInKernel[1])                         = tile(IFFT_X + 5, IFFT_Y + 1 );
            //location<kernel>(ifftGraph.m_r2Comb[1])                               = tile(IFFT_X + 6, IFFT_Y + 1 );
            //location<kernel>(ifftGraph.m_combOutKernel[1])                        = tile(IFFT_X + 7, IFFT_Y + 1 );
    
            ++ifft_graph_insts;
        }
};

class HPGraph: public graph
{
    public:
        input_gmio gmio_in_A[TP_SSR];
        input_gmio gmio_in_B[TP_SSR];
        output_gmio gmio_out[TP_SSR];
      
        HPGraph() {
            dsplib::hadamard::hadamard_graph<TT_DATA_A,
                                             TT_DATA_B,
                                             TP_DIM,
                                             TP_NUM_FRAMES,
                                             TP_HP_SHIFT,
                                             TP_HP_API,
                                             TP_SSR,
                                             TP_RND,
                                             TP_SAT> hpGraph;
         
            runtime<ratio>(*hpGraph.getKernels()) = 0.8;

            for (int i = 0; i < TP_SSR; i++) {
                gmio_in_A[i] = input_gmio::create("gmio_in_vecA_" + std::to_string(hp_graph_insts) + "_" + std::to_string(i), 256, 1000);
                gmio_in_B[i] = input_gmio::create("gmio_in_vecB_" + std::to_string(hp_graph_insts) + "_" + std::to_string(i), 256, 1000);

                connect<>(gmio_in_A[i].out[0], hpGraph.inA[i]);
                connect<>(gmio_in_B[i].out[0], hpGraph.inB[i]);

                gmio_out[i] = output_gmio::create("gmio_out_vec_" + std::to_string(hp_graph_insts) + "_" + std::to_string(i), 256, 1000);
                connect<>(hpGraph.out[i], gmio_out[i].in[0]);
            }

            ++hp_graph_insts;
        }
};

//class CplxConjGraph: public graph
//{
//    private:
//        kernel k_m;
//
//    public:
//        adf::output_gmio gmio_out;
//        adf::input_gmio gmio_in;
//
//        CplxConjGraph() {
//            k_m = adf::kernel::create(cplx_conj_kern);
//            gmio_out = adf::output_gmio::create("gmio_out_cplx_conj_" + std::to_string(cplx_conj_graph_insts), 256, 1000);
//            gmio_in = adf::input_gmio::create("gmio_in_cplx_conj_" + std::to_string(cplx_conj_graph_insts), 256, 1000);
//
//            adf::connect<>(gmio_in.out[0], k_m.in[0]);
//            adf::connect<>(k_m.out[0], gmio_out.in[0]);
//
//            adf::source(k_m) = "cplx_conj.cc";
//
//            adf::runtime<adf::ratio>(k_m) = 0.9;
//
//            ++cplx_conj_graph_insts;
//        }
//};

//class FFTColsGraph: public graph
//{
//    public:
//        input_plio plio_in;
//        output_gmio gmio_out;
//      
//        FFTColsGraph() {
//            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
//                                                         TT_TWIDDLE,
//                                                         TP_POINT_SIZE,
//                                                         TP_FFT,
//                                                         TP_FFT_SHIFT,
//                                                         TP_FFT_CASC_LEN,
//                                                         TP_DYN_PT_SIZE,
//                                                         TP_FFT_WINDOW_VSIZE> fftColsGraph;
//
//            runtime<ratio>(*fftColsGraph.getKernels()) = 0.8;
//
//            // These are port names into and out of the AI Engine
//            std::string data_file_str = std::string(TO_STRING(TT_DATA)) + "/transposed_" + 
//                                        std::to_string(MAT_ROWS) + "x" + 
//                                        std::to_string(MAT_COLS) + 
//                                        "/fft_rows_transposed.txt";
//            plio_in = input_plio::create("cols_plio_in_fft_" + std::to_string(fft_cols_graph_insts), plio_128_bits, data_file_str.c_str());
//            gmio_out = output_gmio::create("cols_gmio_out_fft_" + std::to_string(fft_cols_graph_insts), 128, 1000);
//
//            connect< window<FFT_WINDOW_BUFF_SIZE> > (plio_in.out[0],  fftColsGraph.in[0]);
//            connect< window<FFT_WINDOW_BUFF_SIZE> > (fftColsGraph.out[0], gmio_out.in[0]);
//
//            ++fft_cols_graph_insts;
//        }
//};

//class IFFTRowsGraph: public graph
//{
//    public:
//        input_plio plio_in;
//        output_gmio gmio_out;
//        
//        IFFTRowsGraph() {
//            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
//                                                         TT_TWIDDLE,
//                                                         TP_POINT_SIZE,
//                                                         TP_IFFT,
//                                                         TP_FFT_SHIFT,
//                                                         TP_FFT_CASC_LEN,
//                                                         TP_DYN_PT_SIZE,
//                                                         TP_FFT_WINDOW_VSIZE> ifftRowsGraph;
//           
//            runtime<ratio>(*ifftRowsGraph.getKernels()) = 0.8;
//    
//            // These are port names into and out of the AI Engine
//            std::string data_file_str = std::string(TO_STRING(TT_DATA)) + "/transposed_" + 
//                                        std::to_string(MAT_ROWS) + "x" + 
//                                        std::to_string(MAT_COLS) + 
//                                        "/ifft_cols_transposed.txt";
//            plio_in = input_plio::create("rows_plio_in_ifft_" + std::to_string(ifft_rows_graph_insts), plio_128_bits, data_file_str.c_str());
//            gmio_out = output_gmio::create("rows_gmio_out_ifft_" + std::to_string(ifft_rows_graph_insts), 128, 1000);
//    
//            connect< window<FFT_WINDOW_BUFF_SIZE> > (plio_in.out[0],  ifftRowsGraph.in[0]);
//            connect< window<FFT_WINDOW_BUFF_SIZE> > (ifftRowsGraph.out[0], gmio_out.in[0]);
//    
//            ++ifft_rows_graph_insts;
//        }
//};
//
//
//class PeakGraph: public graph
//{
//    private:
//        kernel cols_k_m;
//        kernel rows_k_m;
//
//    public:
//        adf::output_gmio gmio_out;
//        adf::input_gmio gmio_in;
//
//        PeakGraph() {
//            cols_k_m = adf::kernel::create(cols_peak_kern);
//            rows_k_m = adf::kernel::create(rows_peak_kern);
//            gmio_out = adf::output_gmio::create("gmio_out_peak_" + std::to_string(peak_graph_insts), 256, 1000);
//            gmio_in = adf::input_gmio::create("gmio_in_peak_" + std::to_string(peak_graph_insts), 256, 1000);
//
//            adf::connect<>(gmio_in.out[0], cols_k_m.in[0]);
//            adf::connect<>(cols_k_m.out[0], rows_k_m.in[0]);
//            adf::connect<>(rows_k_m.out[0], gmio_out.in[0]);
//
//            adf::source(cols_k_m) = "peak.cc";
//            adf::source(rows_k_m) = "peak.cc";
//
//            adf::runtime<adf::ratio>(cols_k_m) = 0.9;
//            adf::runtime<adf::ratio>(rows_k_m) = 0.9;
//
//            ++peak_graph_insts;
//        }
//};
