// By: Austin Owens
// Date: 5/23/2024
// Desc: Application to kick off the backprojection engine on the AI cores

#include <unistd.h>
#include <fstream>
#include "sar_backproject.h"
#include <thread>
#include <chrono>

using namespace adf;

int main(int argc, char ** argv) {
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0] << " <xclbin file> <slowtime dataset csv file> <range compressed dataset csv file> <img out csv file> <iteration>" << std::endl;
        return -1;
    }
    std::cout << "\n##### STARTING APPLICATION #####" << std::endl;

    // User args
    char* xclbin_filename = argv[1];
    char* st_dataset_filename = argv[2];
    char* rc_dataset_filename = argv[3];
    char* img_out_filename = argv[4];
    int iter = std::stoi(argv[5]);

    // Instances of the design instantiated across the PL and AIE
    const int INSTANCES = 1;

    std::cout << "\nLoading xclbin, instantiate AIE graph handles, and malloc memory (HOST)..." << std::endl;
    SARBackproject::startTime();
    SARBackproject ifcc(xclbin_filename, st_dataset_filename, rc_dataset_filename, img_out_filename, iter, INSTANCES);
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Init completed (HOST)");

    std::cout << "\nPopulating data buffers (HOST)..." << std::endl;
    SARBackproject::startTime();
    if(ifcc.fetchRadarData() != 0) {
        std::cout << "\nPopulating data buffers failed (HOST)" << std::endl;
        return -1;
    }
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Populating data buffers completed (HOST)");

    std::cout << "\nGenerating target pixels (HOST)..." << std::endl;
    SARBackproject::startTime();
    ifcc.genTargetPixels();
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Generating target pixels  completed (HOST)");

    // Start all AIE kernels
    std::cout << "\nRun AIE graphs (HOST)... " << std::endl;
    SARBackproject::startTime();
    ifcc.runGraphs();
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Run AIE graphs completed (HOST)");

    // Initialize bp xrt::aie::bo buffer arrays 
    xrt::aie::bo buffers_broadcast_data_in[INSTANCES];
    xrt::aie::bo buffers_rc_in[INSTANCES];
    xrt::bo buffers_img_out[INSTANCES];

    // Buffer number for telling AIE kernels how many buffers to operate on in parallel
    int buff_num = 1;

    std::cout << "\nPerform Backprojection (AIE)..." << std::endl;
    buffers_broadcast_data_in[0] = ifcc.m_broadcast_data_buffer;
    buffers_rc_in[0] = ifcc.m_rc_buffer;
    buffers_img_out[0] = ifcc.m_img_buffers[0];
    SARBackproject::startTime();
    ifcc.bp(buffers_broadcast_data_in, buffers_rc_in, 
            buffers_img_out, buff_num);
    SARBackproject::endTime();
    SARBackproject::printTimeDiff("Backprojection completed (AIE)");
    
    if(ifcc.writeImg() != 0) {
        std::cout << "\nWriting image failed!" << std::endl;
        return -1;
    }

    return 0;
};
