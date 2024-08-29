// By: Austin Owens
// Date: 7/1/2024
// Desc: Custom DMA stride controller

#include "dma_hls.h"

////////////////////////////////////////////////////////////
// Top Function of HLS DMA stride controller
////////////////////////////////////////////////////////////
int dma_hls(ap_uint<64>* ddr_mem,
    hls::stream<ap_axiu<128, 0, 0, 0>> &stream_out)
{
    #pragma HLS INTERFACE axis port=stream_out
    int mat_size = 1048576;
    
    #pragma HLS INTERFACE m_axi port=ddr_mem offset=slave bundle=gmem
    #pragma HLS INTERFACE s_axilite port=ddr_mem bundle=control
    #pragma HLS INTERFACE s_axilite port=return bundle=control  
    
    int stride = 1024;
    ap_axiu<128, 0, 0, 0> data_out;
    STRIDE_OUTER:for(int start=0; start < stride; start++) {
        STRIDE_INNER:for(int idx = start; idx < mat_size; idx += stride*2) {
            data_out.data.range(63, 0) = ddr_mem[idx];
            data_out.data.range(127, 64) = ddr_mem[idx+stride];
            data_out.keep = -1;
            stream_out.write(data_out);
        }
    }

    return 0;
}
