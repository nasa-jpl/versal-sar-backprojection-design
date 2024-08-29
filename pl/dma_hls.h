// By: Austin Owens
// Date: 7/3/2024
// Desc: Header for custom DMA controller that handles transpose operations

#pragma once

#include <ap_int.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include "../common.h"

//#if DATA_TYPE == 0
//    using ddr_mem_type = ap_uint<32>[TP_POINT_SIZE][TP_POINT_SIZE];
//#else
//    using ddr_mem_type = ap_uint<64>[TP_POINT_SIZE][TP_POINT_SIZE];
//#endif

////////////////////////////////////////////////////////////
// Top Function of HLS DMA stride controller
////////////////////////////////////////////////////////////
int dma_hls(
    ap_uint<64>* ddr_mem,
    hls::stream<ap_axiu<128, 0, 0, 0>> &stream_out);
