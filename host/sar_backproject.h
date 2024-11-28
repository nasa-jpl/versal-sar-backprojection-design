// By: Austin Owens
// Date: 6/27/2024
// Desc: Header for class to handle SAR backprojection
#ifndef SARBACKPROJECT_H
#define SARBACKPROJECT_H

#include <time.h>
#include <adf.h>
#include <hdf5.h>
#include "xrt/xrt_kernel.h"
#include "xrt/xrt_graph.h"
#include "xrt/xrt_aie.h"
#include "../common.h"

class SARBackproject {
    private:
        // Xclbin and hdf5 filename
        const char* m_xclbin_filename;
        const char* m_hdf5_filename;

        // Number of iterations for AIE kernels to run. The acutal number of
        // iterations of the AIE kernels will be some multiple of m_iter.
        const int m_iter;

        // Instances of the design instantiated across the PL and AIE
        const int m_instances;

        // Device handler and uuid to access acceleration resources
        xrt::device m_device;
        xrt::uuid m_uuid;

        // Define HDF5 data types
        hid_t file;

        // AIE Graphs
        std::vector<xrt::graph> m_bp_graph_hdls;
        //std::vector<xrt::graph> m_fft_graph_hdls;
        //std::vector<xrt::graph> m_ifft_graph_hdls;
        //std::vector<xrt::graph> m_cplx_conj_graph_hdls;
        //std::vector<xrt::graph> m_hp_graph_hdls;

        // PRIVATE FUNCTIONS
        bool readCplxDataset(hid_t file, const std::string& dataset_name, std::vector<TT_DATA>& data);
        bool readFloatDataset(hid_t file, const std::string& dataset_name, std::vector<float>& data);
        

    public:
        // Buffer objects and corresponding arrays for AIE processing
        xrt::aie::bo m_x_ant_pos_buffer;
        float* m_x_ant_pos_array;
        xrt::aie::bo m_y_ant_pos_buffer;
        float* m_y_ant_pos_array;
        xrt::aie::bo m_z_ant_pos_buffer;
        float* m_z_ant_pos_array;
        xrt::aie::bo m_ref_range_buffer;
        float* m_ref_range_array;
        xrt::aie::bo m_xy_px_buffer;
        TT_DATA* m_xy_px_array;
        xrt::aie::bo m_rc_buffer;
        TT_DATA* m_rc_array;
        xrt::aie::bo m_img_buffer;
        TT_DATA* m_img_array;
        //xrt::aie::bo m_range_data_buffer;
        //TT_DATA* m_range_data_array;
        //xrt::aie::bo m_ref_func_buffer;
        //TT_DATA* m_ref_func_array;

    public:
        // Constructor
        SARBackproject(const char* xclbin_filename, const char* h5_filename, int iter, int instances);

        // PUBLIC FUNCTIONS
        static void startTime();
        static void endTime();
        static void resetTimer();
        static void printTimeDiff(const char *msg);
        static void printTotalTime(int curr_iter);
        static void printAvgTime(int iterations);
        static void reshapeMatrix(TT_DATA* input, int rows, int cols, int segment, bool reverse=false);
        static void strideCols(TT_DATA* input, int rows, int cols, int stride, bool reverse=false);
        void runGraphs();
        bool fetchRadarData();
        void bp(xrt::aie::bo* buffers_x_ant_pos_in, xrt::aie::bo* buffers_y_ant_pos_in, 
                xrt::aie::bo* buffers_z_ant_pos_in, xrt::aie::bo* buffers_ref_range_in,
                xrt::aie::bo* buffers_xy_px_in, xrt::aie::bo* buffers_rc_in, 
                xrt::aie::bo* buffers_img_out, int num_of_buffers);
        void fft(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers);
        void ifft(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers);
        void cplxConj(xrt::aie::bo* buffers_in, xrt::aie::bo* buffers_out, int num_of_buffers);
        void elemMatMult(xrt::aie::bo* buffersA_in, xrt::aie::bo* buffersB_in, 
                         xrt::aie::bo* buffers_out, int num_of_buffers);

        // PUBLIC INSTANCE VARIABLES
        // Block size of images
        static constexpr int BLOCK_SIZE_ENTRIES = MAT_ROWS * MAT_COLS;
        static constexpr int BLOCK_SIZE_BYTES = BLOCK_SIZE_ENTRIES * sizeof(TT_DATA);

        // PL dma_hls fft kernel handlers and buffers 
        //std::vector<xrt::kernel> m_dma_hls_fft_kernels;
        //std::vector<xrt::run> m_dma_hls_fft_run_hdls;
        //std::vector<xrt::bo> m_dma_hls_fft_buffers;
        //std::vector<xrt::aie::bo> m_aie_to_pl_buffers;

        // PL dma_hls ifft kernel handlers and buffers 
        //std::vector<xrt::kernel> m_dma_hls_ifft_kernels;
        //std::vector<xrt::run> m_dma_hls_ifft_run_hdls;
        //std::vector<xrt::bo> m_dma_hls_ifft_buffers;

        // Static time metric vars
        static double total_time;
        static double total_avg_time;
        static struct timespec time_start;
        static struct timespec time_end;

};

#endif // SARBACKPROJECT_H
