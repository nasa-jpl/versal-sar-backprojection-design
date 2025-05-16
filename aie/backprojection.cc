// By: Austin Owens
// Date: 10/16/2024
// Desc: Performs backprojection image reconstruction on desired pixel targets across multiple pulses

#include "aie_api/aie_adf.hpp"
#include "aie_api/utils.hpp"
#include "custom_kernels.h"

using namespace adf;

// Supports up to 32 stream connections
void px_arbiter_kern(input_stream<float>* __restrict px_xyz_in,
                     output_pktstream *px_xyz_out) {
    

    // 3-bit pattern to identify the type of packet
    int pkt_type = 0;
    
    // Pixel components per AI reconstruction kernel
    int px_components_per_ai = ((PULSES*RC_SAMPLES)/IMG_SOLVERS)*3;

    for (int switch_idx=0; switch_idx<IMG_SOLVERS; switch_idx++) {
        // Get packet ID for routing from specific index. Packet ID is automatically
        // given at compile time and must be fetched indirectly via an index.
        uint32 pkt_id = getPacketid(px_xyz_out, switch_idx);
        writeHeader(px_xyz_out, pkt_type, pkt_id);

        for (int px_idx=0; px_idx<px_components_per_ai; px_idx++) {
            float px_xyz_val = readincr(px_xyz_in);
        
            //printf("%d: switch_idx: %d | pkt_id: %d | px_xyz_val: %f | tlast: %d\n", 
            //        px_idx, switch_idx, pkt_id, px_xyz_val, px_idx==(16*3)-1);

            writeincr(px_xyz_out, px_xyz_val, px_idx==(px_components_per_ai-1));
        }
    }
}

void slowtime_splicer_kern(input_buffer<float, extents<1>>& __restrict x_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict y_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict z_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict ref_range_in,
                           input_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_in,
                           output_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_out,
                           output_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_out) {
    // Antenna X, Y, and Z position and reference range to scene center
    auto x_ant_pos_in_iter = aie::begin(x_ant_pos_in);
    auto y_ant_pos_in_iter = aie::begin(y_ant_pos_in);
    auto z_ant_pos_in_iter = aie::begin(z_ant_pos_in);
    auto ref_range_in_iter = aie::begin(ref_range_in);

    // Range compressed samples
    auto rc_in_iter = aie::begin_vector<16>(rc_in);

    // Output to backprojection kernel(s)
    auto st_out_iter = aie::begin(slowtime_out);
    auto rc_out_iter = aie::begin_vector<16>(rc_out);
    
    for(unsigned i=0; i < 1; i++) chess_prepare_for_pipelining {
        *st_out_iter++ = *x_ant_pos_in_iter++;
        *st_out_iter++ = *y_ant_pos_in_iter++;
        *st_out_iter++ = *z_ant_pos_in_iter++;
        *st_out_iter++ = *ref_range_in_iter++;
    }

    for(unsigned i=0; i < RC_SAMPLES/16; i++) chess_prepare_for_pipelining {
        *rc_out_iter++ = *rc_in_iter++;
    }
}

ImgReconstruct::ImgReconstruct(int id)
: m_id(id)
//, m_iter(0)
//, m_prev_low_rc((cfloat){0.0f, 0.0f})
//, m_prev_img_val((cfloat){0.0f, 0.0f})
{}

void ImgReconstruct::img_reconstruct_kern(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
                                          input_async_buffer<cfloat, extents<RC_SAMPLES>>& __restrict rc_in,
                                          input_pktstream *px_xyz_in,
                                          //output_pktstream *img_out,
                                          output_async_buffer<cfloat, extents<(PULSES*RC_SAMPLES)/IMG_SOLVERS>>& __restrict img_out,
                                          int rtp_dump_img_in) {

    //printf("%d: rtp_rc_idx_offset_in: %d\n", m_id, rtp_rc_idx_offset_in);

    // Lock range compressed async buffer
    rc_in.acquire();
    //img_out.acquire();
    
    const int SAMPLES = (PULSES*RC_SAMPLES)/IMG_SOLVERS;

    // Used for buffering X, Y, and Z target pixels from input
    // NOTE: Target pixels are stored as ints since thats the datatype associated
    // with pktswitch streams
    alignas(aie::vector_decl_align) float x_trgt_pxls[SAMPLES];
    alignas(aie::vector_decl_align) float y_trgt_pxls[SAMPLES];
    alignas(aie::vector_decl_align) float z_trgt_pxls[SAMPLES];

    // Extract target pixels into buffer for prcessing later in this function
    uint32 header = readincr(px_xyz_in);
    //uint32 id = header & 0x1F;
    //uint32 pkt_type = (header & 0x7000) >> 12;
    //printf("%d IMR_RECON: header: 0x%08x | id: %u | pkt_type: %u\n", m_id, header, id, pkt_type);
    for(int px_idx=0; px_idx < SAMPLES; px_idx++) chess_prepare_for_pipelining {
        int raw_x = readincr(px_xyz_in);
        int raw_y = readincr(px_xyz_in);
        int raw_z = readincr(px_xyz_in);

        x_trgt_pxls[px_idx] = *reinterpret_cast<float*>(&raw_x);
        y_trgt_pxls[px_idx] = *reinterpret_cast<float*>(&raw_y);
        z_trgt_pxls[px_idx] = *reinterpret_cast<float*>(&raw_z);
    }
    
    // Initalize m_img (don't need to do this if keeping track of valid bounds via RTP params)
    //auto zero_init_vec = aie::zeros<cfloat, 16>();
    //for(int i=0; i<RC_SAMPLES/IMG_SOLVERS/16; i++) {
    //    zero_init_vec.store(m_img + i*16);
    //}

    // Initialize radar params
    float ph_corr_coef = (4*PI*MIN_FREQ)/C;

    //// Initialize range and azimuth grid
    int half_size = RC_SAMPLES/2;
    
    // Extract antenna position and range to scene center from slow time data
    auto st_in_vec_iter = aie::begin_vector<ST_ELEMENTS>(slowtime_in);
    auto st_in_vec = *st_in_vec_iter++;
    float x_ant = st_in_vec[0];
    float y_ant = st_in_vec[1];
    float z_ant = st_in_vec[2];
    float r0 = st_in_vec[3];
    
    // Ranged compressed data in time domain (compressed range lines)
    auto rc_in_iter = aie::begin(rc_in);

    // Image output
    auto img_out_iter = aie::begin_vector<16>(img_out);

    // Declare X, Y, and Z target pixel int vectors
    //aie::vector<int,16> x_pxls_int_vec; //= aie::zeros<int32,16>();
    //aie::vector<int,16> y_pxls_int_vec; //= aie::zeros<int32,16>();
    //aie::vector<int,16> z_pxls_int_vec; //= aie::zeros<int32,16>();

    // Declare X, Y, and Z target pixel float vectors
    aie::vector<float,16> x_pxls_vec;
    aie::vector<float,16> y_pxls_vec;
    aie::vector<float,16> z_pxls_vec;
    
    for(int px_seg_idx=0; px_seg_idx < SAMPLES/16; px_seg_idx++) chess_prepare_for_pipelining {

        // Extract X, Y and Z target pixel buffers into their respective pixel vectors
        x_pxls_vec.load(x_trgt_pxls + px_seg_idx*16);
        y_pxls_vec.load(y_trgt_pxls + px_seg_idx*16);
        z_pxls_vec.load(z_trgt_pxls + px_seg_idx*16);

        //uint32 header = readincr(px_xyz_in);
        //uint32 id = header & 0x1F;
        //uint32 pkt_type = (header & 0x7000) >> 12;
        ////printf("%d IMR_RECON: header: 0x%08x | id: %u | pkt_type: %u\n", m_id, header, id, pkt_type);
        //
        //// Extract x, y and z pixel vectors and convert from int32 to float
        //aie::vector<int32,16> x_pxls_int_vec = aie::zeros<int32,16>();
        //aie::vector<int32,16> y_pxls_int_vec = aie::zeros<int32,16>();
        //aie::vector<int32,16> z_pxls_int_vec = aie::zeros<int32,16>();

        //for(int i=0; i<16; i++) {
        //    x_pxls_int_vec.set(readincr(px_xyz_in), i);
        //    y_pxls_int_vec.set(readincr(px_xyz_in), i);
        //    z_pxls_int_vec.set(readincr(px_xyz_in), i);
        //}

        //aie::vector<float,16> x_pxls_vec = x_pxls_int_vec.cast_to<float>();
        //aie::vector<float,16> y_pxls_vec = y_pxls_int_vec.cast_to<float>();
        //aie::vector<float,16> z_pxls_vec = z_pxls_int_vec.cast_to<float>();


        //auto xy_px_vec = *xy_in_iter++;
        //auto x_pxls_vec = aie::real(xy_px_vec);
        //auto y_pxls_vec = aie::imag(xy_px_vec);
        //auto z_pxls_vec = *z_in_iter++;

        printf("%d IMG_RECON: x_pxls_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
                x_pxls_vec.get(0),  x_pxls_vec.get(1),  x_pxls_vec.get(2),  x_pxls_vec.get(3), 
                x_pxls_vec.get(4),  x_pxls_vec.get(5),  x_pxls_vec.get(6),  x_pxls_vec.get(7),
                x_pxls_vec.get(8),  x_pxls_vec.get(9),  x_pxls_vec.get(10), x_pxls_vec.get(11), 
                x_pxls_vec.get(12), x_pxls_vec.get(13), x_pxls_vec.get(14), x_pxls_vec.get(15));
        printf("%d IMG_RECON: y_pxls_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
                y_pxls_vec.get(0),  y_pxls_vec.get(1),  y_pxls_vec.get(2),  y_pxls_vec.get(3), 
                y_pxls_vec.get(4),  y_pxls_vec.get(5),  y_pxls_vec.get(6),  y_pxls_vec.get(7),
                y_pxls_vec.get(8),  y_pxls_vec.get(9),  y_pxls_vec.get(10), y_pxls_vec.get(11), 
                y_pxls_vec.get(12), y_pxls_vec.get(13), y_pxls_vec.get(14), y_pxls_vec.get(15));
        
        /**** CALCULATE DIFFERENTIAL RANGE FOR PIXEL SEGMENTS ****/
        // Calculate this on host to corr pixels to rc...calc boundaries/extremes instead of the whole thing on host
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
        //printf("%d IMG_RECON: differ_range_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //        differ_range_vec.get(0),  differ_range_vec.get(1),  differ_range_vec.get(2),  differ_range_vec.get(3), 
        //        differ_range_vec.get(4),  differ_range_vec.get(5),  differ_range_vec.get(6),  differ_range_vec.get(7),
        //        differ_range_vec.get(8),  differ_range_vec.get(9),  differ_range_vec.get(10), differ_range_vec.get(11), 
        //        differ_range_vec.get(12), differ_range_vec.get(13), differ_range_vec.get(14), differ_range_vec.get(15));

        // Calculate the approximate index for differ_range_vec in the equally spaced range grid
        auto px_idx_acc = aie::mul(INV_RANGE_RES, differ_range_vec);

        // Shift indices to align with proper range compressed samples
        px_idx_acc = aie::add(px_idx_acc, (float)(half_size));

        // Round to nearest whole number
        auto low_idx_float_vec = aie::sub(px_idx_acc, 0.5f).to_vector<float>(0);
        auto low_idx_int_vec = aie::to_fixed<int32>(low_idx_float_vec);
        auto high_idx_int_vec = aie::add(low_idx_int_vec, 1);

        //printf("%d IMG_RECON: px_bounds={%d, %d}\n", m_id, low_px_idx_bound, high_px_idx_bound);
        //printf("%d IMG_RECON: low_idx_int_vec=[%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d]\n", m_id, 
        //        low_idx_int_vec.get(0),  low_idx_int_vec.get(1),  low_idx_int_vec.get(2),  low_idx_int_vec.get(3), 
        //        low_idx_int_vec.get(4),  low_idx_int_vec.get(5),  low_idx_int_vec.get(6),  low_idx_int_vec.get(7),
        //        low_idx_int_vec.get(8),  low_idx_int_vec.get(9),  low_idx_int_vec.get(10), low_idx_int_vec.get(11), 
        //        low_idx_int_vec.get(12), low_idx_int_vec.get(13), low_idx_int_vec.get(14), low_idx_int_vec.get(15));

        // Determine proper boundary of pixel to see if this AI kernel instance even 
        // has the range compressed samples to support the rest of the calculation.
        // Higher values in the pixel vector = lower value pixel indices
        //if ((low_idx_int_vec.get(15) > high_px_idx_bound) || (high_idx_int_vec.get(0) < low_px_idx_bound)) {
        //    //x_pxls_vec = aie::add(x_pxls_vec, range_res*16);
        //    printf("%d: skip: %d > %d || %d < %d\n", m_id, low_idx_int_vec.get(15), high_px_idx_bound, high_idx_int_vec.get(0), low_px_idx_bound);
        //    continue;
        //}

        //**** CALCULATE PHASE CORRECTION FOR IMAGE ****//
        auto ph_corr_angle_vec = aie::mul(ph_corr_coef, differ_range_vec).to_vector<float>(0);
        //printf("%d IMG_RECON: ph_corr_coef = %f\n", m_id, ph_corr_coef);
        //printf("%d IMG_RECON: ph_corr_angle_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //        ph_corr_angle_vec.get(0),  
        //        ph_corr_angle_vec.get(1),  
        //        ph_corr_angle_vec.get(2),  
        //        ph_corr_angle_vec.get(3), 
        //        ph_corr_angle_vec.get(4),  
        //        ph_corr_angle_vec.get(5),  
        //        ph_corr_angle_vec.get(6),  
        //        ph_corr_angle_vec.get(7), 
        //        ph_corr_angle_vec.get(8), 
        //        ph_corr_angle_vec.get(9),  
        //        ph_corr_angle_vec.get(10), 
        //        ph_corr_angle_vec.get(11), 
        //        ph_corr_angle_vec.get(12), 
        //        ph_corr_angle_vec.get(13), 
        //        ph_corr_angle_vec.get(14), 
        //        ph_corr_angle_vec.get(15));

        // Figure out the number of times 2*PI goes into ph_corr_angle_vec. 
        // Floor round to neg infinity by casting to int32, then back to float for later operations
        //auto ph_corr_angle_vec_add_pi = aie::add(ph_corr_angle_vec, PI);
        auto num_pi_wrapped_acc = aie::mul(INV_TWO_PI, ph_corr_angle_vec);
        auto num_pi_wrapped_floor_vec = aie::sub(num_pi_wrapped_acc, 0.5f).to_vector<float>(0);
        auto num_pi_wrapped_int_vec = aie::to_fixed<int32>(num_pi_wrapped_floor_vec); // Round to nearest whole number
        num_pi_wrapped_floor_vec = aie::to_float(num_pi_wrapped_int_vec);
        //printf("%d IMG_RECON: num_pi_wrapped_floor_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //        num_pi_wrapped_floor_vec.get(0),  
        //        num_pi_wrapped_floor_vec.get(1),  
        //        num_pi_wrapped_floor_vec.get(2),  
        //        num_pi_wrapped_floor_vec.get(3), 
        //        num_pi_wrapped_floor_vec.get(4),  
        //        num_pi_wrapped_floor_vec.get(5),  
        //        num_pi_wrapped_floor_vec.get(6),  
        //        num_pi_wrapped_floor_vec.get(7), 
        //        num_pi_wrapped_floor_vec.get(8), 
        //        num_pi_wrapped_floor_vec.get(9),  
        //        num_pi_wrapped_floor_vec.get(10), 
        //        num_pi_wrapped_floor_vec.get(11), 
        //        num_pi_wrapped_floor_vec.get(12), 
        //        num_pi_wrapped_floor_vec.get(13), 
        //        num_pi_wrapped_floor_vec.get(14), 
        //        num_pi_wrapped_floor_vec.get(15));

        // Scale down ph_corr_angle to be within valid domain for sin/cos operation (must be between -PI to PI; modulus doesn't exist)
        auto scale_down_angle_acc = aie::negmul(TWO_PI, num_pi_wrapped_floor_vec); //tmp1
        scale_down_angle_acc = aie::sub(scale_down_angle_acc, PI); //tmp2
        //auto tmp = scale_down_angle_acc.to_vector<float>(0);
        //printf("%d IMG_RECON: REDUCED scale_down_angle_acc=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //        tmp.get(0),  
        //        tmp.get(1),  
        //        tmp.get(2),  
        //        tmp.get(3), 
        //        tmp.get(4),  
        //        tmp.get(5),  
        //        tmp.get(6),  
        //        tmp.get(7), 
        //        tmp.get(8), 
        //        tmp.get(9),  
        //        tmp.get(10), 
        //        tmp.get(11), 
        //        tmp.get(12), 
        //        tmp.get(13), 
        //        tmp.get(14), 
        //        tmp.get(15));
        ph_corr_angle_vec = aie::add(scale_down_angle_acc, ph_corr_angle_vec).to_vector<float>(0); //tmp3
        //printf("%d IMG_RECON: REDUCED ph_corr_angle_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //        ph_corr_angle_vec.get(0),  
        //        ph_corr_angle_vec.get(1),  
        //        ph_corr_angle_vec.get(2),  
        //        ph_corr_angle_vec.get(3), 
        //        ph_corr_angle_vec.get(4),  
        //        ph_corr_angle_vec.get(5),  
        //        ph_corr_angle_vec.get(6),  
        //        ph_corr_angle_vec.get(7), 
        //        ph_corr_angle_vec.get(8), 
        //        ph_corr_angle_vec.get(9),  
        //        ph_corr_angle_vec.get(10), 
        //        ph_corr_angle_vec.get(11), 
        //        ph_corr_angle_vec.get(12), 
        //        ph_corr_angle_vec.get(13), 
        //        ph_corr_angle_vec.get(14), 
        //        ph_corr_angle_vec.get(15));
        //phCorr = exp(1i*4*pi*(data.minF(ii)/c)*dR)
        //float ph_corr_coef = (4*PI*min_freq)/C;

        // Calculate the sin and cos of ph_corr_angle and store as a cfloat (cos in the real part, sin in the imaginary)
        auto ph_corr_vec = aie::sincos_complex(ph_corr_angle_vec);
        ph_corr_vec = aie::neg(ph_corr_vec);
        auto ph_corr_real = aie::real(ph_corr_vec);
        auto ph_corr_imag = aie::imag(ph_corr_vec);
        //printf("%d IMG_RECON: ph_corr=[{%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}, {%f, %f}]\n", m_id, 
        //        ph_corr_real.get(0),  
        //        ph_corr_imag.get(0),  
        //        ph_corr_real.get(1),  
        //        ph_corr_imag.get(1),  
        //        ph_corr_real.get(2),  
        //        ph_corr_imag.get(2),  
        //        ph_corr_real.get(3), 
        //        ph_corr_imag.get(3),  
        //        ph_corr_real.get(4),  
        //        ph_corr_imag.get(4),  
        //        ph_corr_real.get(5),  
        //        ph_corr_imag.get(5),  
        //        ph_corr_real.get(6),  
        //        ph_corr_imag.get(6),  
        //        ph_corr_real.get(7), 
        //        ph_corr_imag.get(7),  
        //        ph_corr_real.get(8), 
        //        ph_corr_imag.get(8),  
        //        ph_corr_real.get(9),  
        //        ph_corr_imag.get(9),  
        //        ph_corr_real.get(10), 
        //        ph_corr_imag.get(10),  
        //        ph_corr_real.get(11), 
        //        ph_corr_imag.get(11),  
        //        ph_corr_real.get(12), 
        //        ph_corr_imag.get(12),  
        //        ph_corr_real.get(13), 
        //        ph_corr_imag.get(13),  
        //        ph_corr_real.get(14), 
        //        ph_corr_imag.get(14),  
        //        ph_corr_real.get(15),
        //        ph_corr_imag.get(15)); 
        
        // Fractional part for interpolation
        low_idx_float_vec = aie::to_float(low_idx_int_vec);
        auto px_delta_idx_vec = aie::sub(px_idx_acc, low_idx_float_vec).to_vector<float>(0);
        
        // TODO: ASSUMING PASSING ALL RC DATA IN
        ////auto shifted_high_idx_int_vec = high_idx_int_vec; //aie::sub(high_idx_int_vec, low_px_idx_bound);
        ////auto shifted_low_idx_int_vec = low_idx_int_vec; //aie::sub(low_idx_int_vec, low_px_idx_bound);
        //auto shifted_high_idx_int_vec = aie::sub(high_idx_int_vec, rtp_rc_idx_offset_in);
        //auto shifted_low_idx_int_vec = aie::sub(low_idx_int_vec, rtp_rc_idx_offset_in);
        //printf("%d IMG_RECON: shifted_low_idx_int_vec=[%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d]\n", m_id, 
        //        shifted_low_idx_int_vec.get(0),  shifted_low_idx_int_vec.get(1),  shifted_low_idx_int_vec.get(2),  shifted_low_idx_int_vec.get(3), 
        //        shifted_low_idx_int_vec.get(4),  shifted_low_idx_int_vec.get(5),  shifted_low_idx_int_vec.get(6),  shifted_low_idx_int_vec.get(7),
        //        shifted_low_idx_int_vec.get(8),  shifted_low_idx_int_vec.get(9),  shifted_low_idx_int_vec.get(10), shifted_low_idx_int_vec.get(11), 
        //        shifted_low_idx_int_vec.get(12), shifted_low_idx_int_vec.get(13), shifted_low_idx_int_vec.get(14), shifted_low_idx_int_vec.get(15));
        cfloat rc_delta;
        for(int px_idx=0; px_idx<16; px_idx++) chess_prepare_for_pipelining {
            //if (low_idx_int_vec.get(i) >= low_px_idx_bound && high_idx_int_vec.get(i) <= high_px_idx_bound) {
            //printf("%d: shifted_low_idx_int_vec.get(%d): %d\n", m_id, i, shifted_low_idx_int_vec.get(i));
            //
            // If we have been in this loop before, write the previous
            // image val out to next kernel in chain
            //if(prev_in_loop) {
            //    writeincr(img_out, m_prev_img_val);
            //    prev_in_loop = false;
            //} else {
            //    // Send low rc boundary to next kernel in chain if not the last
            //    if (m_id != IMG_SOLVERS-1)
            //        writeincr(img_out, rc_in_iter[0]);

            //    // If this is not the first image reconstruction kernel, 
            //    // send all contents from prior kernel to next kernel.
            //    if (m_id != 0) {
            //        bool tlast = false;

            //        // Extract low rc boundary from previous kernel in chain
            //        prev_kern_low_rc_bound = readincr(img_in, tlast);

            //        // Pass values obtained from previous kernel to next kernel in chain 
            //        while(!tlast) {
            //            writeincr(img_out, readincr(img_in, tlast));
            //        }
            //    }
            //}

            //// If this is true, that means we need to use the previous rc segment 
            //// data. This is needed because the interpolation needs to be made between 
            //// two diffrent rc segments. 
            //if (high_idx_int_vec.get(i) == high_px_idx_bound)
            //    rc_delta = prev_kern_low_rc_bound - rc_in_iter[shifted_low_idx_int_vec.get(i)];
            //else
            //    rc_delta = rc_in_iter[shifted_high_idx_int_vec.get(i)] - rc_in_iter[shifted_low_idx_int_vec.get(i)];

            //rc_delta = rc_in_iter[shifted_high_idx_int_vec.get(px_idx)] - rc_in_iter[shifted_low_idx_int_vec.get(px_idx)];
            rc_delta = rc_in_iter[high_idx_int_vec.get(px_idx)] - rc_in_iter[low_idx_int_vec.get(px_idx)];
            auto px_rc_delta = rc_delta*(float)(px_delta_idx_vec.get(px_idx));
            //auto interp = px_rc_delta+rc_in_iter[shifted_low_idx_int_vec.get(px_idx)];
            auto interp = px_rc_delta+rc_in_iter[low_idx_int_vec.get(px_idx)];

            auto img = interp*ph_corr_vec.get(px_idx);
            //auto img = interp;
            //printf("%d: idx=%d\n", m_id, (m_id*(RC_SAMPLES/IMG_SOLVERS))+(px_seg_idx*16 + i));
            
            //printf("%d: img[%d] = %f, %f\n", m_id, (px_seg_idx*16) + i, aie::real(img), aie::imag(img));
            m_img[(px_seg_idx*16) + px_idx] += img;
            //printf("%d: m_img[%d] = %f, %f\n", m_id, (px_seg_idx*16) + i, m_img[(px_seg_idx*16) + i].real, m_img[(px_seg_idx*16) + i].imag);
            
            // Save previous low rc data in case it's needed for the next rc segment
            //m_prev_low_rc = rc_in_iter[shifted_low_idx_int_vec.get(i)];
            
            // Save previous image value so it can be pushed to img_out axi 
            // stream on next iteration. Need to do it this way to properly 
            // handle when to send TLAST.
            //m_prev_img_val = img;

            //m_img[i] = img;
            //prev_in_loop = true;

            //elem_cnt++;
            //}
        }
        //x_pxls_vec = aie::add(x_pxls_vec, range_res*16);
        
        //// Write out TLAST on last image element
        //if (prev_in_loop)
        //    writeincr(img_out, m_prev_img_val, true);

        //y_pxls_vec = aie::sub(y_pxls_vec, az_res);
        //
        //// Start back at beginning of x pixels
        //x_pxls_vec.load(x_pxls_start);
        //prev_in_loop = false;
    }


    //rtp_img_elem_cnt_out = elem_cnt;
    //m_iter++;
    
    // TODO: I THINK THE IMG_OUT_ITER OUTPUT IS READING FALSE DATA BECAUSE THE PRINT IS RIGHT.
    if (rtp_dump_img_in) {

        //auto img_iter = aie::begin(m_img);
        auto img_iter = aie::begin_vector<16>(m_img);
        //for(int i=0; i<SAMPLES; i++) {
        for(int i=0; i<SAMPLES/16; i++) {
            *img_out_iter++ = *img_iter++;
            //cfloat val = *img_iter++;
            //writeincr(img_out, (int)val.real, false);
            //writeincr(img_out, (int)val.imag, i==(SAMPLES-1));

        }

    }

    // Release buffers to show this kernel is finished working on them
    rc_in.release();
    img_out.release();
}


