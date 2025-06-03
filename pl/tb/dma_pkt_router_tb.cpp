// By: Austin Owens
// Date: 5/24/2025
// Desc: Test DMA packet router Programmable Logic (PL) kernel. The DMA packet
// router is an AI packet routing switch that undoes the randomization caused
// by AI pipelining optimization.

#include "../dma_pkt_router.h"
#include <fstream>
#include <iostream>

using namespace std;

// It helps to imagine this testbench from the AI Engine perspective

int main() {
    // Total number of AXIS 128b samples. The divide by 2 comes from the AXIS
    // bus being 128 bits (one cfloat is 64 bits) and the +IMG_SOLVERS comes
    // from each img solver kernel generating a 128b metadata header
    const int AXIS128_SAMPLES = (PULSES*RC_SAMPLES)/2 + IMG_SOLVERS;

    hls::stream<ap_axiu<128, 0, 0, 0>> aie_stream_in;

    // Allocate DDR memory buffer
    ap_uint<64>* ddr_mem = (ap_uint<64>*) malloc(PULSES*RC_SAMPLES*8);
    for (int i = 0; i < PULSES*RC_SAMPLES; i++)
        ddr_mem[i] = 0;
    
    // Loop through AIE data AXI Streams (each iteration is representative of calling another instance of the PL kernel)
    for(int pl_kern=0; pl_kern<AIE_SWITCHES; pl_kern++) {

        // Current working dir inside versal-design-build/build/hw/plsim/dma_pkt_router_testbench/solution1/csim/build
        // Assumes there is only one instance of bpGraph (hence the 0 on aie_to_plio_switch_0_*)
        std::string aie_data_str = "../../../../../aiesim/aiesimulator_output/aie_to_plio_switch_0_" + std::to_string(pl_kern) + ".csv";
        std::ifstream infile(aie_data_str.c_str());
        if (!infile.is_open()) {
            std::cerr << "ERROR: Cannot open input CSV file.\n\n";
            return 1;
        }

        // Read and discard the header line
        std::string header;
        std::getline(infile, header);

        // Parse CSV lines into the stream
        std::string line;
        while (std::getline(infile, line)) {
            // skip non-DATA lines
            if (line.rfind("DATA:", 0) != 0) 
                continue;

            //printf("%s\n", line.c_str());
            
            std::stringstream ss(line);
            std::string token;
            std::vector<std::string> tokens;
            while (std::getline(ss, token, ',')) {
                // trim whitespace
                token.erase(0, token.find_first_not_of(" \t"));
                token.erase(token.find_last_not_of(" \t") + 1);
                tokens.push_back(token);
            }

            // Parse the four 32-bit data fields (D0, D1, D2, D3)
            uint32_t d0, d1, d2, d3;
            d0 = static_cast<uint32_t>(std::stoul(tokens[1]));
            d1 = static_cast<uint32_t>(std::stoul(tokens[2]));
            d2 = static_cast<uint32_t>(std::stoul(tokens[3]));
            d3 = static_cast<uint32_t>(std::stoul(tokens[4]));

            // Parse TLAST (0 or 1)
            int tlast_flag = std::stoi(tokens[5]);

            // Parse TKEEP (which might be -1 or hex). -1 indicates all bytes valid.
            std::string tkeep_str = tokens[6];
            unsigned tkeep_val;
            if (tkeep_str == "-1") {
                // -1 in CSV means all valid bytes (0xFFFF for 128-bit bus)
                tkeep_val = 0xFFFF;
            } else if (tkeep_str.rfind("0x", 0) == 0 || tkeep_str.rfind("0X", 0) == 0) {
                tkeep_val = std::stoul(tkeep_str, nullptr, 16);
            } else if (!tkeep_str.empty()) {
                tkeep_val = std::stoul(tkeep_str);
            } else {
                // Empty TKEEP means no effect, treat as full usage
                tkeep_val = 0xFFFF;
            }

            // We can ignore TIME_NS for functionality  it's only a timestamp
            // double time_ns = std::stod(tokens[7]);  // not used further
            
            // Populate the data128 with the 128 bits of data
            ap_uint<128> data128 = 0;
            data128.range(31,  0)  = d0;
            data128.range(63, 32)  = d1;
            data128.range(95, 64)  = d2;
            data128.range(127,96)  = d3;

            // Define the AXIS data that the PL kernel will ingest
            ap_axiu<128,0,0,0> pkt;
            pkt.data = data128;
            pkt.last = tlast_flag;
            pkt.keep = (ap_uint<16>) (tkeep_val & 0xFFFF);

            // Write data into aie_stream_in
            aie_stream_in.write(pkt);
        }

        // PL kernel invocation
        int ret = dma_pkt_router(aie_stream_in, ddr_mem);

    }

    // Open a file for writing. Current working dir inside 
    // versal-design-build/design/pl/tb/dma_pkt_router_testbench/solution1/csim/build
    FILE *img_fp = fopen("../../../../output_img.csv", "w");
    if (img_fp == NULL) {
        perror("Error opening output_img.csv file");
        return 1;
    }

    float *ddr_mem_f = reinterpret_cast<float*>(ddr_mem);

    fprintf(img_fp, "%.12f%+.12fi", ddr_mem_f[0], ddr_mem_f[1]);
    for(int i=1; i<PULSES*RC_SAMPLES; i++) {
        if (i%RC_SAMPLES == 0) {
            fprintf(img_fp, "\n");
        }
        fprintf(img_fp, ",%.12f%+.12fi", ddr_mem_f[i*2], ddr_mem_f[i*2 + 1]);
    }


    // Print output buffer contents
    printf("Raw AIE AXI Stream data (identical to aiesim CSV content but correctly reordered):\n");
    for (int i = 0; i < PULSES*RC_SAMPLES; i++)
        printf("ddr_mem[%d] = %u, %u\n", i, 
                                         (uint32_t)(ddr_mem[i] & 0xFFFFFFFF),
                                         (uint32_t)(ddr_mem[i]>>32));

    printf("\nSuccessfully wrote output image to versal-design-build/plsim/output_img.csv\n\n");

    free(ddr_mem);

    return 0;
}
