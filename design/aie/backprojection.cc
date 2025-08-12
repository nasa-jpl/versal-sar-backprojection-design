// By: Austin Owens
// Date: 10/16/2024
// Desc: Performs backprojection image reconstruction on desired pixel targets
// across multiple pulses

#include "aie_api/aie_adf.hpp"
#include "aie_api/utils.hpp"
#include "custom_kernels.h"

using namespace adf;

// Supports up to 32 stream connections
void px_demux_kern(input_stream<float>* __restrict px_xyz_in,
                   output_pktstream *px_xyz_out) {
    

    // 3-bit pattern to identify the type of packet
    int pkt_type = 0;
    
    aie::vector<float,4> px_xyz_vec;

    for (int pulse_idx=0; pulse_idx<PULSES; pulse_idx++) chess_prepare_for_pipelining {
        for (int switch_idx=0; switch_idx<IMG_SOLVERS_PER_SWITCH; switch_idx++) chess_prepare_for_pipelining {

            // Get packet ID for routing from specific index. Packet ID is
            // automatically given at compile time and must be fetched
            // indirectly via an index.
            uint32 pkt_id = getPacketid(px_xyz_out, switch_idx);

            writeHeader(px_xyz_out, pkt_type, pkt_id);

            for (int i=0; i<BLOCK_SIZE; i++) chess_prepare_for_pipelining {
                px_xyz_vec = readincr_v<4>(px_xyz_in);
            
                writeincr(px_xyz_out, px_xyz_vec.get(0));
                writeincr(px_xyz_out, px_xyz_vec.get(1));
                writeincr(px_xyz_out, px_xyz_vec.get(2));
                writeincr(px_xyz_out, px_xyz_vec.get(3), i==(BLOCK_SIZE-1));
            }
        }
    }
}

void data_broadcast_kern(input_stream<float>* __restrict slowtime_in,
                         input_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_in,
                         output_stream<float>* __restrict slowtime_out,
                         output_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_out) {

    // Input and output iterators for range compressed samples
    auto rc_in_iter = aie::begin_vector<16>(rc_in);
    auto rc_out_iter = aie::begin_vector<16>(rc_out);
    
    // Write out slowtime data
    writeincr(slowtime_out, readincr_v<4>(slowtime_in));
    
    // Write out RC data
    for(unsigned i=0; i < RC_SAMPLES/16; i++) chess_prepare_for_pipelining {
        *rc_out_iter++ = *rc_in_iter++;
    }
}

ImgReconstruct::ImgReconstruct(int id)
: m_id(id)
{}

void ImgReconstruct::img_reconstruct_kern(input_buffer<float, extents<BC_ELEMENTS>>& __restrict slowtime_in,
                                          input_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_in,
                                          input_pktstream *px_xyz_in,
                                          output_pktstream *img_out,
                                          int rtp_dump_img_in) {

    const int SAMPLES = (AZ_SAMPLES*RC_SAMPLES)/IMG_SOLVERS;

    // Initialize radar params
    float ph_corr_coef = (4*PI*MIN_FREQ)/C;

    // Initialize range and azimuth grid
    int half_size = RC_SAMPLES/2;
    
    // Extract antenna position and range to scene center from slow time data
    auto st_in_vec_iter = aie::begin_vector<BC_ELEMENTS>(slowtime_in);
    auto st_in_vec = *st_in_vec_iter++;
    float x_ant = st_in_vec[0];
    float y_ant = st_in_vec[1];
    float z_ant = st_in_vec[2];
    float r0 = st_in_vec[3];
    
    // Ranged compressed data in time domain (compressed range lines)
    auto rc_in_iter = aie::begin(rc_in);

    // Declare X, Y, and Z target pixel float vectors
    aie::vector<float,16> x_pxls_vec;
    aie::vector<float,16> y_pxls_vec;
    aie::vector<float,16> z_pxls_vec;

    int raw_x, raw_y, raw_z;
    for(int px_seg_idx=0; px_seg_idx < SAMPLES/16; px_seg_idx++) chess_prepare_for_pipelining {

        uint32 header = readincr(px_xyz_in);
        for(int i=0; i<16; i++) chess_prepare_for_pipelining {
            raw_x = readincr(px_xyz_in);
            raw_y = readincr(px_xyz_in);
            raw_z = readincr(px_xyz_in);
            readincr(px_xyz_in); // zero pad
            x_pxls_vec.set(*reinterpret_cast<float*>(&raw_x), i);
            y_pxls_vec.set(*reinterpret_cast<float*>(&raw_y), i);
            z_pxls_vec.set(*reinterpret_cast<float*>(&raw_z), i);
        }

        auto x_vec = aie::sub(x_ant, x_pxls_vec);
        auto x_diff_sq_acc = aie::mul_square(x_vec);

        auto y_vec = aie::sub(y_ant, y_pxls_vec);
        auto y_diff_sq_acc = aie::mul_square(y_vec);

        auto z_vec = aie::sub(z_ant, z_pxls_vec);
        auto z_diff_sq_acc = aie::mul_square(z_vec);

        auto xy_add_acc = aie::add(x_diff_sq_acc, y_diff_sq_acc.to_vector<float>(0));
        auto xyz_add_acc = aie::add(xy_add_acc, z_diff_sq_acc.to_vector<float>(0));

        auto xyz_sqrt_acc = aie::sqrt(xyz_add_acc);

        auto differ_range_vec = aie::sub(xyz_sqrt_acc, r0).to_vector<float>(0);

        // Calculate the approximate index for differ_range_vec in the equally
        // spaced range grid
        auto px_idx_acc = aie::mul(INV_RANGE_RES, differ_range_vec);

        // Shift indices to align with proper range compressed samples
        px_idx_acc = aie::add(px_idx_acc, (float)(half_size));

        // Round to nearest whole number
        auto low_idx_float_vec = aie::sub(px_idx_acc, 0.5f).to_vector<float>(0);
        auto low_idx_int_vec = aie::to_fixed<int32>(low_idx_float_vec);
        auto high_idx_int_vec = aie::add(low_idx_int_vec, 1);

        //**** CALCULATE PHASE CORRECTION FOR IMAGE ****//
        auto ph_corr_angle_vec = aie::mul(ph_corr_coef, differ_range_vec).to_vector<float>(0);

        // Figure out the number of times 2*PI goes into ph_corr_angle_vec.
        // Floor round to neg infinity by casting to int32, then back to float
        // for later operations
        auto num_pi_wrapped_acc = aie::mul(INV_TWO_PI, ph_corr_angle_vec);
        auto num_pi_wrapped_floor_vec = aie::sub(num_pi_wrapped_acc, 0.5f).to_vector<float>(0);
        auto num_pi_wrapped_int_vec = aie::to_fixed<int32>(num_pi_wrapped_floor_vec); // Round to nearest whole number
        num_pi_wrapped_floor_vec = aie::to_float(num_pi_wrapped_int_vec);

        // Scale down ph_corr_angle to be within valid domain for sin/cos
        // operation (must be between -PI to PI; modulus doesn't exist)
        auto scale_down_angle_acc = aie::negmul(TWO_PI, num_pi_wrapped_floor_vec);
        scale_down_angle_acc = aie::sub(scale_down_angle_acc, PI);
        ph_corr_angle_vec = aie::add(scale_down_angle_acc, ph_corr_angle_vec).to_vector<float>(0);

        // Calculate the sin and cos of ph_corr_angle and store as a cfloat
        // (cos in the real part, sin in the imaginary)
        auto ph_corr_vec = aie::sincos_complex(ph_corr_angle_vec);
        ph_corr_vec = aie::neg(ph_corr_vec);
        auto ph_corr_real = aie::real(ph_corr_vec);
        auto ph_corr_imag = aie::imag(ph_corr_vec);

        // Fractional part for interpolation
        low_idx_float_vec = aie::to_float(low_idx_int_vec);
        auto px_delta_idx_vec = aie::sub(px_idx_acc, low_idx_float_vec).to_vector<float>(0);
        
        cfloat rc_delta;
        for(int px_idx=0; px_idx<16; px_idx++) chess_prepare_for_pipelining {
            rc_delta = rc_in_iter[high_idx_int_vec.get(px_idx)] - rc_in_iter[low_idx_int_vec.get(px_idx)];
            auto px_rc_delta = rc_delta*(float)(px_delta_idx_vec.get(px_idx));
            auto interp = px_rc_delta+rc_in_iter[low_idx_int_vec.get(px_idx)];

            auto img = interp*ph_corr_vec.get(px_idx);
            m_img[(px_seg_idx*16) + px_idx] += img;
        }
    }

    if (rtp_dump_img_in) {
        
        // Create iterator for m_img buffer
        auto img_iter = aie::begin(m_img);

        // Get first element from m_img buffer
        cfloat img_elem = *img_iter++;

        // Construct and place uint32 header on img_out stream
        writeHeader(img_out, 0, getPacketid(img_out, 0));

        // Add additional metadata indicating kernel instantiation number
        writeincr(img_out, m_id, false);
    
        // Add two uint32 padding since bus width is 128 bit 
        writeincr(img_out, 0, false);
        writeincr(img_out, 0, false);

        for(int i=0; i<SAMPLES-1; i++) chess_prepare_for_pipelining {
            writeincr(img_out, img_elem.real, false);
            writeincr(img_out, img_elem.imag, false);
            img_elem = *img_iter++;
        }
        writeincr(img_out, img_elem.real, false);
        writeincr(img_out, img_elem.imag, true);
    }
}
