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
        // Xclbin filename
        const char* m_xclbin_filename;
        const char* m_st_dataset_filename;
        const char* m_rc_dataset_filename;
        const char* m_img_out_filename;

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

    public:
        // Buffer objects and corresponding arrays for AIE processing
        xrt::aie::bo m_broadcast_data_buffer;
        float* m_broadcast_data_array;
        xrt::aie::bo m_rc_buffer;
        cfloat* m_rc_array;

    public:
        // Constructor
        SARBackproject(const char* xclbin_filename, 
                       const char* st_dataset_filename, 
                       const char* rc_dataset_filename, 
                       const char* img_out_filename, 
                       int iter, 
                       int instances);

        // PUBLIC FUNCTIONS
        static void startTime();
        static void endTime();
        static void resetTimer();
        static void printTimeDiff(const char *msg);
        static void printTotalTime(int curr_iter);
        static void printAvgTime(int iterations);
        static void unwrap(double* angles);
        bool writeImg();
        bool fetchRadarData();
        void genTargetPixels();
        void runGraphs();
        void bp(xrt::aie::bo* buffers_broadcast_data_in, xrt::aie::bo* buffers_rc_in, 
                xrt::bo* buffers_img_out, int num_of_buffers);

        // PL dma_stride kernel handlers and buffers 
        std::vector<xrt::kernel> m_dma_stride_kernels;
        std::vector<xrt::run> m_dma_stride_run_hdls;
        std::vector<xrt::bo> m_xyz_px_buffers;
        std::vector<uint64_t*> m_xyz_px_arrays;

        // PL dma_pkt_router kernel handlers and buffers 
        std::vector<xrt::kernel> m_dma_pkt_router_kernels;
        std::vector<xrt::run> m_dma_pkt_router_run_hdls;
        std::vector<xrt::bo> m_img_buffers;
        std::vector<cfloat*> m_img_arrays;

        // Static time metric vars
        static double total_time;
        static double total_avg_time;
        static struct timespec time_start;
        static struct timespec time_end;

};

#endif // SARBACKPROJECT_H
