// By: Austin Owens
// Date: 7/3/2024
// Desc: Header for DMA packet router controller

#pragma once

#include <ap_int.h>
#include <hls_stream.h>
#include <ap_axi_sdata.h>
#include "../common.h"

////////////////////////////////////////////////////////////
// Top Function of DMA packet router controller
////////////////////////////////////////////////////////////
int dma_pkt_router(hls::stream<ap_axiu<128, 0, 0, 0>> &aie_stream_in,
                   ap_uint<64>* ddr_mem);
