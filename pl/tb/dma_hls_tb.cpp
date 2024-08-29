#include "../dma_hls.h"
#include <iostream>
#include "../../common.h"

using namespace std;

// It helps to imagine this testbench from the AI Engine perspective

#if DATA_TYPE == 0 // cint16 {1, 1}
    #define GOLDEN_DATA 0x0001000100010001
#else // cfloat {1.5, 1.5}
    #define GOLDEN_DATA 0x3fc000003fc00000
#endif

int main() {
    hls::stream<ap_axiu<128, 0, 0, 0>> aie_stream_out;
    //hls::stream<ap_axiu<128, 0, 0, 0>> stream_trans_fft_plio_out;
    //ap_uint<64> ddr_mem[1024*1024];
    ap_uint<64>* ddr_mem = (ap_uint<64>*) malloc(1024*1024*8);
    //ap_uint<64> ddr_mem[1024*1024] = {0};
    for (int i=0; i<1024*1024; i++) {
        ddr_mem[i] = 0xFFFFFFFFFFFFFFFF - (1024*1024) + i;
    }
    //cout << ddr_mem[0] << " " << ddr_mem[1024] << endl;
    //cout << ddr_mem[1] << " " << ddr_mem[1025] << endl;

    //hls::stream<ap_axiu<128, 0, 0, 0>> stream_trans_ifft_plio_in;
    //hls::stream<ap_axiu<128, 0, 0, 0>> stream_trans_ifft_plio_out;
    #if DATA_TYPE == 0 // cint16
        //int rows = 1024/4;
        //int cols = 1024/4;
        int mat_size = (1024*1024)/4;
    #else // cfloat
        //int rows = 1024/2;
        //int cols = 1024/2;
        int mat_size = (1024*1024)/2;
    #endif

    // PL KERNEL INVOCATION
    dma_hls(ddr_mem, aie_stream_out);
           
    ap_axiu<128, 0, 0, 0> transpose_out;
    cout << "AIE FFT TRANSPOSED DATA INPUT:" << endl;
    for (int i = 0; i < mat_size; ++i) {

        int row = i >> 9; // Divide by 512
        int col = i & 511; // Modulus of 512

        // AIE INPUT FROM PL
        aie_stream_out.read(transpose_out);

        if ((i <= 2000) || (i >= mat_size-2000))
        {
            cout << i << " transpose_out (" << row << "," << col << "): ";
            #if DATA_TYPE == 0 // cint16
                ap_uint<32> val1 = transpose_out.data.range(31, 0);
                ap_uint<32> val2 = transpose_out.data.range(63, 32);
                ap_uint<32> val3 = transpose_out.data.range(95, 64);
                ap_uint<32> val4 = transpose_out.data.range(127, 96);
                cout << val1 << " " << val2 << " " << val3 << " " << val4 << endl;
            #else // cfloat
                ap_uint<64> val1 = transpose_out.data.range(63, 0);
                ap_uint<64> val2 = transpose_out.data.range(127, 64);
                cout << val1 << " " << val2 << endl;
            #endif

        }

    }

    return 0;

}
