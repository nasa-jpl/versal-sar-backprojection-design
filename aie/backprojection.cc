// By: Austin Owens
// Date: 10/16/2024
// Desc: Performs backprojection image reconstruction on desired pixel targets across multiple pulses

#include "aie_api/aie_adf.hpp"
#include "aie_api/utils.hpp"
#include "custom_kernels.h"

using namespace adf;

void slowtime_splicer_kern(input_buffer<float, extents<1>>& __restrict x_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict y_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict z_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict ref_range_in,
                           output_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_out) {
    // Antenna X, Y, and Z position and reference range to scene center
    auto x_ant_pos_in_iter = aie::begin(x_ant_pos_in);
    auto y_ant_pos_in_iter = aie::begin(y_ant_pos_in);
    auto z_ant_pos_in_iter = aie::begin(z_ant_pos_in);
    auto ref_range_in_iter = aie::begin(ref_range_in);

    // Output to backprojection kernel(s)
    auto st_out_iter = aie::begin(slowtime_out);
    
    *st_out_iter++ = *x_ant_pos_in_iter++;
    *st_out_iter++ = *y_ant_pos_in_iter++;
    *st_out_iter++ = *z_ant_pos_in_iter++;
    *st_out_iter++ = *ref_range_in_iter++;
}

// Only needed to satisfy a port connection to img_reconstruct_kern 
void dummy_kern(output_stream<TT_DATA>* __restrict img_out) {}

ImgReconstruct::ImgReconstruct(int id)
: m_id(id)
, m_iter(0)
, m_prev_low_rc((cfloat){0.0f, 0.0f})
, m_prev_img_val((cfloat){0.0f, 0.0f})
{}

void ImgReconstruct::img_reconstruct_kern(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
                                          input_async_buffer<TT_DATA, extents<TP_POINT_SIZE/4>>& __restrict rc_in,
                                          input_stream<TT_DATA>* __restrict img_in, 
                                          output_stream<TT_DATA>* __restrict img_out,
                                          int32 &rtp_img_elem_cnt_out) {
    
    // Setup header for packet switcher
    //uint32 switch_id = getPacketid(img_out, 0);
    //printf("switch_id: %d\n", switch_id);
    //writeHeader(img_out, 0, switch_id);

    // Lock range compressed async buffer
    rc_in.acquire();

    // Initalize RTP values
    //rtp_header_idx_out = -1;

    // Initalize m_img (don't need to do this if keeping track of valid bounds via RTP params)
    //auto zero_init_vec = aie::zeros<TT_DATA, 16>();
    //for(int i=0; i<TP_POINT_SIZE/IMG_SOLVERS/16; i++) {
    //    zero_init_vec.store(m_img + i*16);
    //}

    // Initialize radar params
    float min_freq = 9288080400.0;
    float ph_corr_coef = (4*PI*min_freq)/C;

    float range_freq_step = 1471301.6;
    float range_width = C/(2.0*range_freq_step);
    float range_res = range_width/TP_POINT_SIZE;
    float inv_range_res = 1.0/range_res;

    float total_az = 0.01726837;
    float az_res = C/(2.0*total_az*min_freq);
    float az_width = AZ_POINT_SIZE/az_res;
    //float az_width = 108.4105;

    // Initialize range and azimuth grid
    int half_size = TP_POINT_SIZE/2;
    int seg_size = TP_POINT_SIZE/IMG_SOLVERS;

    int low_px_idx_bound = (seg_size*3)-(m_id*seg_size);
    int high_px_idx_bound = (seg_size*4)-(m_id*seg_size);

    float start_range_px = -range_width/2;
    float end_range_px = range_width/2;
    float start_az_px = az_width/2.0;
    float end_az_px = -az_width/2.0;

    //printf("start_range_px = %f | end_range_px = %f | start_az_px = %f | end_az_px = %f\n", start_range_px, end_range_px, start_az_px, end_az_px);

    aie::vector<float,16> y_pxls_vec = aie::broadcast<float,16>(start_az_px);
    //aie::vector<float,16> y_pxls_vec = aie::broadcast<float,16>(0.0f);
    aie::vector<float,16> x_pxls_vec;
    for(int i=0; i<16; i++) {
        x_pxls_vec.set(start_range_px + i*range_res, i);
    }

    // Save x_pxls_vec since we will be needing to loop back to start for all y axis points
    alignas(aie::vector_decl_align) float x_pxls_start[16];
    x_pxls_vec.store(x_pxls_start);

    //printf("%d IMG_RECON: x_pxls_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
    //        x_pxls_vec.get(0),  x_pxls_vec.get(1),  x_pxls_vec.get(2),  x_pxls_vec.get(3), 
    //        x_pxls_vec.get(4),  x_pxls_vec.get(5),  x_pxls_vec.get(6),  x_pxls_vec.get(7),
    //        x_pxls_vec.get(8),  x_pxls_vec.get(9),  x_pxls_vec.get(10), x_pxls_vec.get(11), 
    //        x_pxls_vec.get(12), x_pxls_vec.get(13), x_pxls_vec.get(14), x_pxls_vec.get(15));

    //cfloat xy_px
    //for(int i=0; i<16; i++) {
    //    start_px_bound
    //}
    

    // -(seg_size*2) -> -4096
    // -(seg_size)   -> -2048
    // (seg_size)    ->  2048
    // (seg_size*2)  ->  4096

    

    
    ////Pixel bounds based on segmented range compressed data
    //float start_px_bound = ((range_width/IMG_SOLVERS)*m_id) - range_width/2;
    //float end_px_bound = ((range_width/IMG_SOLVERS)*(m_id+1)) - range_width/2;
    //
    //m_seg(m_id, m_range_grid.range_width);

    //// Pixel X start and end position on grid
    //float x_px_st = -5.0;
    //float x_px_en = 30.0;

    //// Pixel Y start and end position on grid
    //float y_px_st = -5.0;
    //float y_px_en = 5.0;
    //
    //// Pixel X & Y grid length (32x32)
    //int x_px_grid_len = 32;
    //int y_px_grid_len = 32;
    //
    //// Pixel X & Y step sizes
    //float dx = (x_px_en - x_px_st) / (x_px_grid_len - 1);
    //float dy = (y_px_en - y_px_st) / (y_px_grid_len - 1);
    
    // Extract antenna position and range to scene center from slow time data
    auto st_in_vec_iter = aie::begin_vector<ST_ELEMENTS>(slowtime_in);
    
    auto st_in_vec = *st_in_vec_iter++;

    float x_ant = st_in_vec[0];
    float y_ant = st_in_vec[1];
    float z_ant = st_in_vec[2];
    float r0 = st_in_vec[3];

    // Ranged compressed data in time domain (compressed range lines)
    auto rc_in_iter = aie::begin(rc_in);

    // X and Y target pixels
    //auto xy_in_iter = aie::begin_vector<16>(xy_px_in);

    // Image output
    //auto img_out_iter = aie::begin_vector<16>(img_out);
    
    // Precalculate z_diff_sq
    auto z_diff_sq = z_ant*z_ant;
    
    //bool first_valid_px_idx = false;
    bool prev_in_loop = false;
    int elem_cnt = 0;
    TT_DATA prev_kern_low_rc_bound;
    for(int y_idx=0; y_idx < AZ_POINT_SIZE; y_idx++) chess_prepare_for_pipelining {
        for(int x_idx=0; x_idx < TP_POINT_SIZE/16; x_idx++) chess_prepare_for_pipelining {
            
            /**** CALCULATE DIFFERENTIAL RANGE FOR PIXEL SEGMENTS ****/
            // Calculate this on host to corr pixels to rc...calc boundaries/extremes instead of the whole thing on host
            auto x_vec = aie::sub(x_ant, x_pxls_vec);
            auto x_diff_sq_acc = aie::mul_square(x_vec);

            auto y_vec = aie::sub(y_ant, y_pxls_vec);
            auto y_diff_sq_acc = aie::mul_square(y_vec);

            auto xy_add_vec = aie::add(x_diff_sq_acc, z_diff_sq).to_vector<float>(0);
            auto xyz_add_acc = aie::add(y_diff_sq_acc, xy_add_vec);

            auto xyz_sqrt_acc = aie::sqrt(xyz_add_acc);

            auto differ_range_vec = aie::sub(xyz_sqrt_acc, r0).to_vector<float>(0);

            // Calculate the approximate index for differ_range_vec in the equally spaced range grid
            auto px_idx_acc = aie::mul(inv_range_res, differ_range_vec);

            // Shift indices to align with proper range compressed samples
            px_idx_acc = aie::add(px_idx_acc, (float)(half_size));

            // Round to nearest whole number
            auto low_idx_float_vec = aie::sub(px_idx_acc, 0.5f).to_vector<float>(0);
            auto low_idx_int_vec = aie::to_fixed<int32>(low_idx_float_vec);
            auto high_idx_int_vec = aie::add(low_idx_int_vec, 1);

            // Determine proper boundary of pixel to see if this AI kernel instance even 
            // has the range compressed samples to support the rest of the calculation.
            // Higher values in the pixel vector = lower value pixel indices
            if ((low_idx_int_vec.get(15) > high_px_idx_bound) || (high_idx_int_vec.get(0) < low_px_idx_bound)) {
                x_pxls_vec = aie::add(x_pxls_vec, range_res*16);
                continue;
            }

            //**** CALCULATE PHASE CORRECTION FOR IMAGE ****//
            auto ph_corr_angle_vec = aie::mul(ph_corr_coef, differ_range_vec).to_vector<float>(0);

            // Floor round to neg infinity by casting to int32, then back to float for later operations
            auto num_pi_wrapped_acc = aie::mul(INV_TWO_PI, ph_corr_angle_vec);
            auto num_pi_wrapped_floor_vec = aie::sub(num_pi_wrapped_acc, 0.5f).to_vector<float>(0);
            auto num_pi_wrapped_int_vec = aie::to_fixed<int32>(num_pi_wrapped_floor_vec);
            num_pi_wrapped_floor_vec = aie::to_float(num_pi_wrapped_int_vec);
            
            // Scale down ph_corr_angle to be within valid domain for sin/cos operation (must be between -PI to PI; modulus doesn't exist)
            auto scale_down_angle_acc = aie::negmul(TWO_PI, num_pi_wrapped_floor_vec);
            scale_down_angle_acc = aie::sub(scale_down_angle_acc, PI);
            ph_corr_angle_vec = aie::add(scale_down_angle_acc, ph_corr_angle_vec).to_vector<float>(0);

            // Calculate the sin and cos of ph_corr_angle and store as a cfloat (cos in the real part, sin in the imaginary)
            auto ph_corr_vec = aie::sincos_complex(ph_corr_angle_vec);
            ph_corr_vec = aie::neg(ph_corr_vec);
            auto ph_corr_real = aie::real(ph_corr_vec);
            auto ph_corr_imag = aie::imag(ph_corr_vec);
            
            // Fractional part for interpolation
            low_idx_float_vec = aie::to_float(low_idx_int_vec);
            auto px_delta_idx_vec = aie::sub(px_idx_acc, low_idx_float_vec).to_vector<float>(0);

            auto shifted_high_idx_int_vec = aie::sub(high_idx_int_vec, low_px_idx_bound);
            auto shifted_low_idx_int_vec = aie::sub(low_idx_int_vec, low_px_idx_bound);
            cfloat rc_delta;
            for(int i=0; i<16; i++) chess_prepare_for_pipelining {
                if (low_idx_int_vec.get(i) >= low_px_idx_bound && high_idx_int_vec.get(i) <= high_px_idx_bound) {
                    //printf("%d: shifted_low_idx_int_vec.get(%d): %d\n", m_id, i, shifted_low_idx_int_vec.get(i));
                    //if (!first_valid_px_idx) {
                    //    rtp_valid_low_bound_out = x_idx*16 + i;
                    //    first_valid_px_idx = true;
                    //}
                    
                    // If we have been in this loop before, write the previous
                    // image val out to next kernel in chain
                    if(prev_in_loop) {
                        writeincr(img_out, m_prev_img_val);
                        prev_in_loop = false;
                    } else {
                        // Send low rc boundary to next kernel in chain if not the last
                        if (m_id != IMG_SOLVERS-1)
                            writeincr(img_out, rc_in_iter[0]);

                        // If this is not the first image reconstruction kernel, 
                        // send all contents from prior kernel to next kernel.
                        if (m_id != 0) {
                            bool tlast = false;

                            // Extract low rc boundary from previous kernel in chain
                            prev_kern_low_rc_bound = readincr(img_in, tlast);

                            // Pass values obtained from previous kernel to next kernel in chain 
                            while(!tlast) {
                                writeincr(img_out, readincr(img_in, tlast));
                            }
                        }
                    }

                    // If this is true, that means we need to use the previous rc segment 
                    // data. This is needed because the interpolation needs to be made between 
                    // two diffrent rc segments. 
                    if (high_idx_int_vec.get(i) == high_px_idx_bound)
                        rc_delta = prev_kern_low_rc_bound - rc_in_iter[shifted_low_idx_int_vec.get(i)];
                    else
                        rc_delta = rc_in_iter[shifted_high_idx_int_vec.get(i)] - rc_in_iter[shifted_low_idx_int_vec.get(i)];


                    auto px_rc_delta = rc_delta*(float)(px_delta_idx_vec.get(i));
                    auto interp = px_rc_delta+rc_in_iter[shifted_low_idx_int_vec.get(i)];

                    auto img = interp*ph_corr_vec.get(i);
                    //auto img = interp;
                    
                    // Save previous low rc data in case it's needed for the next rc segment
                    //m_prev_low_rc = rc_in_iter[shifted_low_idx_int_vec.get(i)];
                    
                    // Save previous image value so it can be pushed to img_out axi 
                    // stream on next iteration. Need to do it this way to properly 
                    // handle when to send TLAST.
                    m_prev_img_val = img;

                    //m_img[i] = img;
                    prev_in_loop = true;

                    elem_cnt++;
                }
            }
            x_pxls_vec = aie::add(x_pxls_vec, range_res*16);
        }
        
        // Write out TLAST on last image element
        if (prev_in_loop)
            writeincr(img_out, m_prev_img_val, true);

        y_pxls_vec = aie::sub(y_pxls_vec, az_res);
        
        // Start back at beginning of x pixels
        x_pxls_vec.load(x_pxls_start);
        prev_in_loop = false;
    }


    rtp_img_elem_cnt_out = elem_cnt;
    //m_iter++;

    // Release buffer to show this kernel is finished working on them
    rc_in.release();
}


