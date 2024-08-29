// By: Austin Owens
// Date: 5/23/2024
// Desc: Application to kick off Image FFT Cross Correlation across the system

#include <adf.h>
#include <unistd.h>
#include <fstream>
#include "xrt/xrt_kernel.h"
#include "xrt/xrt_graph.h"
#include "xrt/xrt_aie.h"

using namespace adf;

void ref_func(int32* din, int32* dout, int size) {
	for(int i=0; i<size; i++) {
		dout[i] = din[i] + 5;
	}
}

const int BLOCK_SIZE_BYTES = 1024;
const int BLOCK_SIZE_INTS = BLOCK_SIZE_BYTES / sizeof(int32);

int main(int argc, char ** argv) {
    if (argc < 3) {
        std::cerr << "Usage: " << argv[0] << " <xclbin file> <data file>" << std::endl;
        return -1;
    }

    char* xclbin_filename = argv[1];
    char* data_filename = argv[2];

	// Open xclbin
	auto device = xrt::device(0);

    // Load xclbin
	auto uuid = device.load_xclbin(xclbin_filename);

    // Instantiate Buffer Objects
    int mem_group = 0;

	auto din_buffer = xrt::aie::bo(device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, mem_group);
	int* din_array= din_buffer.map<int*>();

	auto dout_buffer = xrt::aie::bo(device, BLOCK_SIZE_BYTES, xrt::bo::flags::normal, mem_group);
	int* dout_array= dout_buffer.map<int*>();

    int32* dout_ref = (int32*) malloc(BLOCK_SIZE_BYTES);
    std::cout << "GMIO::malloc completed" << std::endl;

    // Read data from file into din_array
    std::ifstream data_file(data_filename);
    if (!data_file.is_open()) {
        std::cerr << "Error opening data file: " << data_filename << std::endl;
        return -1;
    }

    for (int i = 0; i < BLOCK_SIZE_INTS; i++) {
        int value;
        if (!(data_file >> value)) {
            std::cerr << "Error reading data file: insufficient data" << std::endl;
            data_file.close();
            return -1;
        }
        din_array[i] = value;
    }
    data_file.close();

    // Invoke AI Engine asynchronously
    int offset = 0;
	auto ghdl=xrt::graph(device, uuid, "gr");
	din_buffer.async("gr.gmio_in", XCL_BO_SYNC_BO_GMIO_TO_AIE, BLOCK_SIZE_BYTES, offset);
    ghdl.run(1);
	auto dout_buffer_run=dout_buffer.async("gr.gmio_out", XCL_BO_SYNC_BO_AIE_TO_GMIO, BLOCK_SIZE_BYTES, offset);

    // Wait for gmio_out to complete
    dout_buffer_run.wait();

    ref_func(din_array, dout_ref, BLOCK_SIZE_INTS);

    int error=0;
    for(int i=0; i<BLOCK_SIZE_INTS; i++){
		if(dout_array[i] != dout_ref[i]){
			std::cout << "ERROR:dout[" << i << "]=" << dout_array[i] << ",gold=" << dout_ref[i] << std::endl;
			error++;
		}
    }
    
    std::cout << "GMIO transactions finished" << std::endl;

	ghdl.end();
    if(error==0){
		std::cout << "TEST PASSED!" << std::endl;
    }else{
		std::cout << "ERROR!" << std::endl;
    }

    return error;
};
