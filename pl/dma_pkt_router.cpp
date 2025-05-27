// By: Austin Owens
// Date: 5/24/2025
// Desc: AI packet routing switch that undoes the randomization caused by AI
// pipelining optimization

#include "dma_pkt_router.h"

////////////////////////////////////////////////////////////
// Top Function of DMA packet router controller
////////////////////////////////////////////////////////////
int dma_pkt_router(hls::stream<ap_axiu<128, 0, 0, 0>> &aie_stream_in,
                   ap_uint<64>* ddr_mem) {
    #pragma HLS INTERFACE axis port=aie_stream_in
    #pragma HLS INTERFACE m_axi port=ddr_mem offset=slave bundle=gmem depth=PULSES*RC_SAMPLES
    #pragma HLS INTERFACE s_axilite port=ddr_mem bundle=control
    #pragma HLS INTERFACE s_axilite port=return  bundle=control
    
    const int SAMPLES_PER_KERN = (PULSES*RC_SAMPLES)/IMG_SOLVERS;

    // Metadata declarations
    ap_axiu<128, 0, 0, 0> metadata;
    ap_uint<32> header;
    ap_uint<5> pkt_id;
    ap_uint<3> pkt_type;
    ap_uint<32> instance_id;

    // Image data declaration
    ap_axiu<128, 0, 0, 0> img_data;
    
    // Loop through every image solver kernel
    IMG_KERNEL_LOOP:for(int kern=0; kern<IMG_SOLVERS; kern++) {

        // Read in first 128 bits of packet containing metadata
        metadata = aie_stream_in.read();
        
        // First 32 bits is the packet switch header (See AMD's UG1079 docs):
        //
        // Packet Bit Fields:
        //  +-------+----------------------------+
        //  | Bits  | Field                      |
        //  +-------+----------------------------+
        //  |   4-0 | Packet ID                  |
        //  |  11-5 | 7'b0000000                 |
        //  | 14-12 | Packet Type                |
        //  |    15 | 1'b0                       |
        //  | 20-16 | Source Row                 |
        //  | 27-21 | Source Column              |
        //  | 30-28 | 3'b000                     |
        //  |    31 | Odd parity of bits[30:0]   |
        //  +-------+----------------------------+
        header = metadata.data.range(31, 0);
        pkt_id = header & 0x1F;
        pkt_type = (header & 0x7000) >> 12;
        
        // Next 32 bits is instance ID of kernel that generated the data
        instance_id = metadata.data.range(63, 32);

        // Bits 64-128 are padding and ignored. In the future, this could be used
        // for additional metadata
        //metadata.data.range(127, 64);
        
        // DDR Offset
        int ddr_offset = instance_id*SAMPLES_PER_KERN;
        
        // Loop through number of samples in stream. Each 64 bit is 1 cfloat sample
        // (i.e. 32 bit float real + 32 bit float imaginary)
        IMG_DATA_LOOP:for(int idx=0; idx<SAMPLES_PER_KERN; idx+=2) {
            img_data = aie_stream_in.read();
            ddr_mem[ddr_offset + idx] = img_data.data.range(63, 0);
            ddr_mem[ddr_offset + (idx+1)] = img_data.data.range(127, 64);
        }
    }

    return 0;
}

