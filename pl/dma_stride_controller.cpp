// By: Austin Owens
// Date: 6/4/2025
// Desc: DMA stride controller that pushes data into the AIE engine in a
// specific order.

#include "dma_stride_controller.h"

////////////////////////////////////////////////////////////
// Top Function of DMA stride controller
////////////////////////////////////////////////////////////
int dma_stride_controller(ap_uint<64>* ddr_mem, 
                          hls::stream<ap_axiu<128, 0, 0, 0>> pl_stream_out[AIE_SWITCHES]) {
    #pragma HLS INTERFACE axis port=pl_stream_out
    #pragma HLS INTERFACE m_axi port=ddr_mem offset=slave bundle=gmem depth=PULSES*AZ_SAMPLES*RC_SAMPLES*2
    #pragma HLS INTERFACE s_axilite port=ddr_mem bundle=control
    #pragma HLS INTERFACE s_axilite port=return  bundle=control

    //const int SAMPLES_PER_KERN = (PULSES*RC_SAMPLES)/IMG_SOLVERS;

    // Metadata declarations
    //ap_axiu<128, 0, 0, 0> metadata;
    //ap_uint<32> header;
    //ap_uint<5> pkt_id;
    //ap_uint<3> pkt_type;
    //ap_uint<32> instance_id;

    // Image data declaration
    //ap_axiu<128, 0, 0, 0> img_data;

    //int block_size = 16;
    //int stride = (PULSES*RC_SAMPLES*2)/IMG_SOLVERS_PER_SWITCH; //80
    //ap_axiu<128, 0, 0, 0> trgt_px_data;
    //
    //int idx = 0;
    //STRIDE_OUTER:for(int base=0; base < block_size*PULSES*2; base+=block_size*2) {
    //    printf("PL: base=%d\n", base);
    //    STRIDE_MIDDLE:for(int start=base; start < stride*IMG_SOLVERS_PER_SWITCH; start+=stride) {
    //        printf("PL: base=%d | start=%d\n", base, start);
    //        idx = start;
    //        STRIDE_INNER:for(int i=start; i < start+block_size; i++) {
    //            trgt_px_data.data.range(63, 0) = ddr_mem[idx++];
    //            trgt_px_data.data.range(127, 64) = ddr_mem[idx++];
    //            float x, y;
    //            uint32_t raw_x = (uint32_t) trgt_px_data.data.range(31,  0);
    //            uint32_t raw_y = (uint32_t) trgt_px_data.data.range(63, 32);
    //            std::memcpy(&x, &raw_x, sizeof(float));
    //            std::memcpy(&y, &raw_y, sizeof(float));
    //            printf("PL: base=%d | start=%d | i=%d | idx=%d | x=%f | y=%f\n", base, start, i, idx, x, y);
    //            trgt_px_data.keep = -1;
    //            trgt_px_data.last = 0;
    //            pl_stream_out.write(trgt_px_data);
    //        }
    //    }
    //}
    
    // Number of target pixels in a single pulse
    const int px_128b = AZ_SAMPLES * RC_SAMPLES;

    // Number of target pixels destined for each switch per pulse
    const int px_per_sw_128b = px_128b/AIE_SWITCHES;

    // Number of target pixels destined for each AIE kernel per pulse
    const int px_per_kern_128b = px_128b/IMG_SOLVERS;
    
    ap_axiu<128, 0, 0, 0> trgt_px_data;
    int idx_64b;
    
    // Iterate over every pulse
    PULSE_LOOP:for(int pulse_idx = 0; pulse_idx < PULSES; pulse_idx++) {

        // Iterate over each switch
        SWITCH_LOOP:for(int sw_num = 0; sw_num < AIE_SWITCHES; sw_num++) {
            // Starting pixel index on each pulse
            int start_px_idx_per_pulse = sw_num*px_per_sw_128b + pulse_idx*px_128b;

            // Strides BLOCK_SIZE amount of pixels
            STRIDE_OUTER:for(int base_128b = start_px_idx_per_pulse; base_128b < px_per_kern_128b + start_px_idx_per_pulse; base_128b += BLOCK_SIZE) {

                // Strides the number of pixels that will be sent to an AIE kernel
                STRIDE_MIDDLE:for(int start_128b = base_128b; start_128b < px_per_sw_128b + start_px_idx_per_pulse; start_128b += px_per_kern_128b) {

                    // Extract BLOCK_SIZE amount of pixels from DDR and place them onto an AXI Stream
                    #pragma HLS PIPELINE II=1
                    STRIDE_INNER:for(int idx_128b = start_128b; idx_128b < start_128b + BLOCK_SIZE; idx_128b++) {

                        // Because idx_128b is a 128 bit sample index representing
                        // a pixel, but ddr_mem is 64bit, we need to multiply by 2
                        idx_64b = idx_128b*2;
                        
                        // Fetch data from DDR
                        trgt_px_data.data.range(63, 0) = ddr_mem[idx_64b];
                        trgt_px_data.data.range(127, 64) = ddr_mem[idx_64b + 1];
                        trgt_px_data.keep = -1;

                        // TLAST is asserted once every BLOCK_SIZE pixels are let through
                        trgt_px_data.last = (idx_128b == (start_128b + BLOCK_SIZE - 1)) ? 1 : 0;
                        
                        // Write out to specific AXI Stream
                        pl_stream_out[sw_num].write(trgt_px_data);
                    }
                }
            }
        }
    }

    //const int block_size = 16;
    //const int px_128b = AZ_SAMPLES * RC_SAMPLES;
    //const int px_per_kern_128b = px_128b/IMG_SOLVERS;
    //const int base_shift_px = block_size * AZ_SAMPLES;

    //ap_axiu<128, 0, 0, 0> trgt_px_data;
    //int idx_64b;

    //STRIDE_OUTER:for(int base_128b = 0; base_128b < base_shift_px; base_128b += block_size) {
    //    STRIDE_MIDDLE:for(int start_128b = base_128b; start_128b < px_128b; start_128b += px_per_kern_128b) {
    //        #pragma HLS PIPELINE II=1
    //        STRIDE_INNER:for(int idx_128b = start_128b; idx_128b < start_128b + block_size; idx_128b++) {
    //            // idx_128b is the 128 bit sample index. To index ddr_mem
    //            // (which is 64bit), mulitply by 2
    //            idx_64b = idx_128b*2;

    //            trgt_px_data.data.range(63, 0) = ddr_mem[idx_64b];
    //            trgt_px_data.data.range(127, 64) = ddr_mem[idx_64b + 1];
    //            trgt_px_data.keep = -1;

    //            // TLAST is asserted once block_size 128b samples are let through
    //            trgt_px_data.last = (idx_128b == (start_128b + block_size - 1)) ? 1 : 0;

    //            pl_stream_out.write(trgt_px_data);
    //        }
    //    }
    //}

    return 0;
}

