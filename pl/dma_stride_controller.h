// By: Austin Owens
// Date: 6/4/2025
// Desc: Header for DMA stride controller. This controller feeds data into AIE.

#pragma once

#include <ap_int.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include "../common.h"

////////////////////////////////////////////////////////////
// Top Function of DMA packet router controller
////////////////////////////////////////////////////////////
int dma_stride_controller(ap_uint<64>* ddr_mem,
                          hls::stream<ap_axiu<128, 0, 0, 0>> pl_stream_out[AIE_SWITCHES]);
