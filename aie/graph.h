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

extern uint8_t fft_rows_graph_insts;
extern uint8_t fft_cols_graph_insts;
extern uint8_t ifft_rows_graph_insts;
extern uint8_t ifft_cols_graph_insts;
extern uint8_t hp_graph_insts;
extern uint8_t cplx_conj_graph_insts;
extern uint8_t peak_graph_insts;

class FFTRowsGraph: public graph
{
    public:
        input_gmio gmio_in;
        output_gmio gmio_out;
      
        FFTRowsGraph() {
            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
                                                         TT_TWIDDLE,
                                                         TP_POINT_SIZE,
                                                         TP_FFT,
                                                         TP_FFT_SHIFT,
                                                         TP_FFT_CASC_LEN,
                                                         TP_DYN_PT_SIZE,
                                                         TP_FFT_WINDOW_VSIZE> fftRowsGraph;

            runtime<ratio>(*fftRowsGraph.getKernels()) = 0.8;

            // These are port names in to and out of the AI Engine
            gmio_in = input_gmio::create("rows_gmio_in_fft_" + std::to_string(fft_rows_graph_insts), 256, 1000);
            gmio_out = output_gmio::create("rows_gmio_out_fft_" + std::to_string(fft_rows_graph_insts), 256, 1000);

            connect< window<FFT_WINDOW_BUFF_SIZE> > (gmio_in.out[0], fftRowsGraph.in[0]);
            connect< window<FFT_WINDOW_BUFF_SIZE> > (fftRowsGraph.out[0], gmio_out.in[0]);

            ++fft_rows_graph_insts;
        }
};

class FFTColsGraph: public graph
{
    public:
        input_plio plio_in;
        output_gmio gmio_out;
      
        FFTColsGraph() {
            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
                                                         TT_TWIDDLE,
                                                         TP_POINT_SIZE,
                                                         TP_FFT,
                                                         TP_FFT_SHIFT,
                                                         TP_FFT_CASC_LEN,
                                                         TP_DYN_PT_SIZE,
                                                         TP_FFT_WINDOW_VSIZE> fftColsGraph;

            runtime<ratio>(*fftColsGraph.getKernels()) = 0.8;

            // These are port names into and out of the AI Engine
            std::string data_file_str = std::string(TO_STRING(TT_DATA)) + "/transposed_" + 
                                        std::to_string(MAT_ROWS) + "x" + 
                                        std::to_string(MAT_COLS) + 
                                        "/fft_rows_transposed.txt";
            plio_in = input_plio::create("cols_plio_in_fft_" + std::to_string(fft_cols_graph_insts), plio_128_bits, data_file_str.c_str());
            gmio_out = output_gmio::create("cols_gmio_out_fft_" + std::to_string(fft_cols_graph_insts), 128, 1000);

            connect< window<FFT_WINDOW_BUFF_SIZE> > (plio_in.out[0],  fftColsGraph.in[0]);
            connect< window<FFT_WINDOW_BUFF_SIZE> > (fftColsGraph.out[0], gmio_out.in[0]);

            ++fft_cols_graph_insts;
        }
};

class IFFTColsGraph: public graph
{
    public:
        input_gmio gmio_in;
        output_gmio gmio_out;
        
        IFFTColsGraph() {
            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
                                                         TT_TWIDDLE,
                                                         TP_POINT_SIZE,
                                                         TP_IFFT,
                                                         TP_FFT_SHIFT,
                                                         TP_FFT_CASC_LEN,
                                                         TP_DYN_PT_SIZE,
                                                         TP_FFT_WINDOW_VSIZE> ifftColsGraph;

            runtime<ratio>(*ifftColsGraph.getKernels()) = 0.8;
    
            // These are port names in to and out of the AI Engine
            gmio_in = input_gmio::create("cols_gmio_in_ifft_" + std::to_string(ifft_cols_graph_insts), 256, 1000);
            gmio_out = output_gmio::create("cols_gmio_out_ifft_" + std::to_string(ifft_cols_graph_insts), 256, 1000);
    
            connect< window<FFT_WINDOW_BUFF_SIZE> > (gmio_in.out[0], ifftColsGraph.in[0]);
            connect< window<FFT_WINDOW_BUFF_SIZE> > (ifftColsGraph.out[0], gmio_out.in[0]);
    
            ++ifft_cols_graph_insts;
        }
};

class IFFTRowsGraph: public graph
{
    public:
        input_plio plio_in;
        output_gmio gmio_out;
        
        IFFTRowsGraph() {
            dsplib::fft::dit_1ch::fft_ifft_dit_1ch_graph<TT_DATA,
                                                         TT_TWIDDLE,
                                                         TP_POINT_SIZE,
                                                         TP_IFFT,
                                                         TP_FFT_SHIFT,
                                                         TP_FFT_CASC_LEN,
                                                         TP_DYN_PT_SIZE,
                                                         TP_FFT_WINDOW_VSIZE> ifftRowsGraph;
           
            runtime<ratio>(*ifftRowsGraph.getKernels()) = 0.8;
    
            // These are port names into and out of the AI Engine
            std::string data_file_str = std::string(TO_STRING(TT_DATA)) + "/transposed_" + 
                                        std::to_string(MAT_ROWS) + "x" + 
                                        std::to_string(MAT_COLS) + 
                                        "/ifft_cols_transposed.txt";
            plio_in = input_plio::create("rows_plio_in_ifft_" + std::to_string(ifft_rows_graph_insts), plio_128_bits, data_file_str.c_str());
            gmio_out = output_gmio::create("rows_gmio_out_ifft_" + std::to_string(ifft_rows_graph_insts), 128, 1000);
    
            connect< window<FFT_WINDOW_BUFF_SIZE> > (plio_in.out[0],  ifftRowsGraph.in[0]);
            connect< window<FFT_WINDOW_BUFF_SIZE> > (ifftRowsGraph.out[0], gmio_out.in[0]);
    
            ++ifft_rows_graph_insts;
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
                                             TP_API,
                                             TP_SSR,
                                             TP_RND,
                                             TP_SAT> hpGraph;
         
        runtime<ratio>(*hpGraph.getKernels()) = 0.8;

        for (int i = 0; i < TP_SSR; i++) {
            gmio_in_A[i] = input_gmio::create("gmio_in_vecA_" + std::to_string(hp_graph_insts) + "_" + std::to_string(i), 256, 1000);
            gmio_in_B[i] = input_gmio::create("gmio_in_vecB_" + std::to_string(hp_graph_insts) + "_" + std::to_string(i), 256, 1000);

            connect<>(gmio_in_A[i].out[0],  hpGraph.inA[i]);
            connect<>(gmio_in_B[i].out[0],  hpGraph.inB[i]);

            gmio_out[i] = output_gmio::create("gmio_out_vec_" + std::to_string(hp_graph_insts) + "_" + std::to_string(i), 256, 1000);
            connect<>(hpGraph.out[i], gmio_out[i].in[0]);
        }

        ++hp_graph_insts;
    }
};

class CplxConjGraph: public graph
{
    private:
        kernel k_m;

    public:
        adf::output_gmio gmio_out;
        adf::input_gmio gmio_in;

        CplxConjGraph() {
            k_m = adf::kernel::create(cplx_conj_kern);
            gmio_out = adf::output_gmio::create("gmio_out_cplx_conj_" + std::to_string(cplx_conj_graph_insts), 256, 1000);
            gmio_in = adf::input_gmio::create("gmio_in_cplx_conj_" + std::to_string(cplx_conj_graph_insts), 256, 1000);

            adf::connect<>(gmio_in.out[0], k_m.in[0]);
            adf::connect<>(k_m.out[0], gmio_out.in[0]);

            adf::source(k_m) = "cplx_conj.cc";

            adf::runtime<adf::ratio>(k_m) = 0.9;

            ++cplx_conj_graph_insts;
        }
};

class PeakGraph: public graph
{
    private:
        kernel cols_k_m;
        kernel rows_k_m;

    public:
        adf::output_gmio gmio_out;
        adf::input_gmio gmio_in;

        PeakGraph() {
            cols_k_m = adf::kernel::create(cols_peak_kern);
            rows_k_m = adf::kernel::create(rows_peak_kern);
            gmio_out = adf::output_gmio::create("gmio_out_peak_" + std::to_string(peak_graph_insts), 256, 1000);
            gmio_in = adf::input_gmio::create("gmio_in_peak_" + std::to_string(peak_graph_insts), 256, 1000);

            adf::connect<>(gmio_in.out[0], cols_k_m.in[0]);
            adf::connect<>(cols_k_m.out[0], rows_k_m.in[0]);
            adf::connect<>(rows_k_m.out[0], gmio_out.in[0]);

            adf::source(cols_k_m) = "peak.cc";
            adf::source(rows_k_m) = "peak.cc";

            adf::runtime<adf::ratio>(cols_k_m) = 0.9;
            adf::runtime<adf::ratio>(rows_k_m) = 0.9;

            ++peak_graph_insts;
        }
};
