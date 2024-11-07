// By: Austin Owens
// Date: 10/16/2024
// Desc: Performs phase correction compressed range lines for target

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

Backprojection::Backprojection(int id)
: m_id(id)
, m_iter(0)
{
    // Initialize m_img to all zeros
    auto img_iter = aie::begin_vector<16>(m_img);
    auto zero_init_vec = aie::zeros<TT_DATA, 16>();
    for(unsigned i=0; i<2048/16; i++) {
        *img_iter++ = zero_init_vec;
    }
}

void Backprojection::init_radar_params(float range_freq_step, int range_samples, float min_freq)
{
    m_radar_params.range_freq_step = range_freq_step;
    m_radar_params.range_samples = range_samples;
    m_radar_params.min_freq = min_freq;
    m_radar_params.ph_corr_coef = (4*PI*min_freq)/C;
    
    // Using 2*PI instead of 4*pi to account for 2*PI unwrapping that happens later in algorithm
    //m_radar_params.ph_corr_coef = (2*PI*min_freq)/C; 
}

void Backprojection::init_range_grid()
{
    m_range_grid.range_width = C/(2.0*m_radar_params.range_freq_step);
    m_range_grid.range_res = m_range_grid.range_width/m_radar_params.range_samples;
    m_range_grid.inv_range_res = 1.0/m_range_grid.range_res;
    m_range_grid.seg_offset = (float)((3-m_id)*2048);
    
    //int sample_num = -2048;
    //for (int i=0; i<BP_SOLVERS; i++) {
    //    int id = i-(BP_SOLVERS/2);
    //    if (id < 0)
    //        id+=BP_SOLVERS;

    //    if(m_id == id) {
    //        m_range_grid.samp_start_size = sample_num;
    //        break;
    //    }

    //    sample_num+=2048;
    //}
}

void Backprojection::print_float(const char *str, aie::vector<float,16> data) {
    printf("%d: %s=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, str,
           data.get(0),  data.get(1),  data.get(2),  data.get(3), 
           data.get(4),  data.get(5),  data.get(6),  data.get(7),
           data.get(8),  data.get(9),  data.get(10), data.get(11), 
           data.get(12), data.get(13), data.get(14), data.get(15));

}

void Backprojection::print_int32(const char *str, aie::vector<int32,16> data) {
    printf("%d: %s=%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", m_id, str,
           data.get(0),  data.get(1),  data.get(2),  data.get(3), 
           data.get(4),  data.get(5),  data.get(6),  data.get(7),
           data.get(8),  data.get(9),  data.get(10), data.get(11), 
           data.get(12), data.get(13), data.get(14), data.get(15));

}

//void Backprojection::init_pixel_grid(float x_st, float x_en, int x_len, float y_st, float y_en, int y_len)
//{
//    m_px_grid.x_st = x_st;
//    m_px_grid.x_en = x_en;
//    m_px_grid.x_len = x_len;
//    m_px_grid.x_res = (x_en - x_st) / (x_len - 1);
//    m_px_grid.y_st = y_st;
//    m_px_grid.y_en = y_en;
//    m_px_grid.y_len = y_len;
//    m_px_grid.y_res = (y_en - y_st) / (y_len - 1);
//}
//
//void Backprojection::init_pixel_segment()
//{
//    m_seg.range_st = ((m_range_grid.range_width/BP_SOLVERS)*m_id) - m_range_grid.range_width/2;
//    m_seg.range_en = ((m_range_grid.range_width/BP_SOLVERS)*(m_id+1)) - m_range_grid.range_width/2;
//
//    m_seg.st_x_px_bound = false;
//    m_seg.en_x_px_bound = false;
//    m_seg.all_x_px_bound = false;
//    if (m_px_grid.x_st >= m_seg.range_st && m_px_grid.x_st < m_seg.range_en) {
//        m_seg.st_x_px_bound = true;
//    }
//    if (m_px_grid.x_en >= m_seg.range_st && m_px_grid.x_en < m_seg.range_en) {
//        m_seg.en_x_px_bound = true;
//    }
//    if (m_px_grid.x_st <= m_seg.range_st && m_px_grid.x_en > m_seg.range_en) {
//        m_seg.all_x_px_bound = true;
//    }
//}

void Backprojection::backprojection_kern(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
                                         input_circular_buffer<TT_DATA, extents<2048>>& __restrict rc_in,
                                         input_buffer<TT_DATA, extents<2048>>& __restrict xy_px_in,
                                         output_async_buffer<TT_DATA, extents<2048>>& __restrict img_out) {
    // Lock img_out buffer
    img_out.acquire();

    // Extract antenna position and range to scene center from slow time data
    auto st_in_vec_iter = aie::begin_vector<ST_ELEMENTS>(slowtime_in);
    auto st_in_vec = *st_in_vec_iter++;

    float x_ant = st_in_vec[0];
    float y_ant = st_in_vec[1];
    float z_ant = st_in_vec[2];
    float r0 = st_in_vec[3];

    // Ranged compressed data in time domain (compressed range lines)
    auto rc_in_iter = aie::begin_random_circular(rc_in);

    // X and Y target pixels
    auto xy_in_iter = aie::begin_vector<16>(xy_px_in);

    // Image output
    auto img_out_iter = aie::begin_vector<16>(img_out);
    
    // Initialize radar parameters
    init_radar_params(1471301.6, 8192, 9288080400.0);

    // Initialize range grid based on range compressed data
    init_range_grid();
    
    // Initialize target pixel grid
    //init_pixel_grid(-5, 30, 32, -5, 5, 32);

    // Initialize pixel segment for this specific kernel
    //init_pixel_segment();

    //Constants
    //const float dF = 1471301.6;
    //int fft_seg = (TP_POINT_SIZE/BP_SOLVERS)*(m_id+1);
    //int fft_seg = (8192/BP_SOLVERS)*(m_id);

    // Calculate max width scene size of image (in range direction)
    //float max_px_width = c/(2*dF);

    //// Range step sizes
    //float dr = max_px_width/8192;
    //float inv_dr = 1.0/dr;
    //m_range_grid(dF, 8192);

    // Pixel bounds based on segmented range compressed data
    //float start_px_bound = ((max_px_width/BP_SOLVERS)*m_id) - max_px_width/2;
    //float end_px_bound = ((max_px_width/BP_SOLVERS)*(m_id+1)) - max_px_width/2;
    
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

    // Image memory
    //alignas(aie::vector_decl_align) TT_DATA img[m_px_grid.x_len*m_px_grid.y_len];
    
    // Precalculate z_diff_sq
    auto z_diff_sq = z_ant*z_ant;

    //for(int xy_idx=0; xy_idx < 2048/16; xy_idx++) chess_prepare_for_pipelining {
    for(int xy_idx=0; xy_idx < 2; xy_idx++) chess_prepare_for_pipelining {
        printf("\n%d: ITERATION: %d\n", m_id, xy_idx);
        
        auto xy_px_vec = *xy_in_iter++;
        auto x_pxls_vec = aie::real(xy_px_vec);
        auto y_pxls_vec = aie::imag(xy_px_vec);
        
        //aie::print(x_pxls, true, "x_pxls=");
        //printf("%d: x_pxls=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
        //        x_pxls.get(0),  x_pxls.get(1),  x_pxls.get(2),  x_pxls.get(3), 
        //        x_pxls.get(4),  x_pxls.get(5),  x_pxls.get(6),  x_pxls.get(7),
        //        x_pxls.get(8),  x_pxls.get(9),  x_pxls.get(10), x_pxls.get(11), 
        //        x_pxls.get(12), x_pxls.get(13), x_pxls.get(14), x_pxls.get(15));

        ////aie::print(y_pxls, true, "y_pxls=");
        //printf("%d: y_pxls=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
        //        y_pxls.get(0),  y_pxls.get(1),  y_pxls.get(2),  y_pxls.get(3), 
        //        y_pxls.get(4),  y_pxls.get(5),  y_pxls.get(6),  y_pxls.get(7),
        //        y_pxls.get(8),  y_pxls.get(9),  y_pxls.get(10), y_pxls.get(11), 
        //        y_pxls.get(12), y_pxls.get(13), y_pxls.get(14), y_pxls.get(15));

        // Calculate differential range for pixel segments
        auto x_vec = aie::sub(x_ant, x_pxls_vec);
        auto x_diff_sq_acc = aie::mul_square(x_vec);

        auto y_vec = aie::sub(y_ant, y_pxls_vec);
        auto y_diff_sq_acc = aie::mul_square(y_vec);

        auto xy_add_vec = aie::add(x_diff_sq_acc, z_diff_sq).to_vector<float>(0);
        auto xyz_add_acc = aie::add(y_diff_sq_acc, xy_add_vec);

        auto xyz_sqrt_acc = aie::sqrt(xyz_add_acc);

        auto differ_range_vec = aie::sub(xyz_sqrt_acc, r0).to_vector<float>(0);
        print_float("differ_range_vec", differ_range_vec);
        //printf("%d: differ_range_vec=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
        //        differ_range_vec.get(0),  differ_range_vec.get(1),  differ_range_vec.get(2),  differ_range_vec.get(3), 
        //        differ_range_vec.get(4),  differ_range_vec.get(5),  differ_range_vec.get(6),  differ_range_vec.get(7),
        //        differ_range_vec.get(8),  differ_range_vec.get(9),  differ_range_vec.get(10), differ_range_vec.get(11), 
        //        differ_range_vec.get(12), differ_range_vec.get(13), differ_range_vec.get(14), differ_range_vec.get(15));

        // Calculate phase correction for image
        auto ph_corr_angle_vec = aie::mul(m_radar_params.ph_corr_coef, differ_range_vec).to_vector<float>(0);
        //print_float("ph_corr_angle_vec", ph_corr_angle_vec);
        //aie::print(ph_corr_angle_vec, true, "ph_corr_angle_vec (before scale)=");
        //printf("%d: ph_corr_angle_vec (before scale)=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //        ph_corr_angle_vec.get(0),  ph_corr_angle_vec.get(1),  ph_corr_angle_vec.get(2),  ph_corr_angle_vec.get(3), 
        //        ph_corr_angle_vec.get(4),  ph_corr_angle_vec.get(5),  ph_corr_angle_vec.get(6),  ph_corr_angle_vec.get(7),
        //        ph_corr_angle_vec.get(8),  ph_corr_angle_vec.get(9),  ph_corr_angle_vec.get(10), ph_corr_angle_vec.get(11), 
        //        ph_corr_angle_vec.get(12), ph_corr_angle_vec.get(13), ph_corr_angle_vec.get(14), ph_corr_angle_vec.get(15));

        // Floor round to neg infinity by casting to int32, then back to float for later operations
        auto num_pi_wrapped_acc = aie::mul(INV_TWO_PI, ph_corr_angle_vec);
        auto num_pi_wrapped_floor_vec = aie::sub(num_pi_wrapped_acc, 0.5f).to_vector<float>(0);
        auto num_pi_wrapped_int_vec = aie::to_fixed<int32>(num_pi_wrapped_floor_vec);
        num_pi_wrapped_floor_vec = aie::to_float(num_pi_wrapped_int_vec);
        //printf("%d: num_pi_wrapped_floor_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //        num_pi_wrapped_floor_vec.get(0),  num_pi_wrapped_floor_vec.get(1),  num_pi_wrapped_floor_vec.get(2),  num_pi_wrapped_floor_vec.get(3), 
        //        num_pi_wrapped_floor_vec.get(4),  num_pi_wrapped_floor_vec.get(5),  num_pi_wrapped_floor_vec.get(6),  num_pi_wrapped_floor_vec.get(7),
        //        num_pi_wrapped_floor_vec.get(8),  num_pi_wrapped_floor_vec.get(9),  num_pi_wrapped_floor_vec.get(10), num_pi_wrapped_floor_vec.get(11), 
        //        num_pi_wrapped_floor_vec.get(12), num_pi_wrapped_floor_vec.get(13), num_pi_wrapped_floor_vec.get(14), num_pi_wrapped_floor_vec.get(15));
        
        // Scale down ph_corr_angle to be within valid domain for sin/cos operation (must be between -PI to PI; modulus doesn't exist)
        auto scale_down_angle_acc = aie::negmul(TWO_PI, num_pi_wrapped_floor_vec);
        scale_down_angle_acc = aie::sub(scale_down_angle_acc, PI);
        ph_corr_angle_vec = aie::add(scale_down_angle_acc, ph_corr_angle_vec).to_vector<float>(0);
        //printf("%d: ph_corr_angle (after scale)=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //        ph_corr_angle_vec.get(0),  ph_corr_angle_vec.get(1),  ph_corr_angle_vec.get(2),  ph_corr_angle_vec.get(3), 
        //        ph_corr_angle_vec.get(4),  ph_corr_angle_vec.get(5),  ph_corr_angle_vec.get(6),  ph_corr_angle_vec.get(7),
        //        ph_corr_angle_vec.get(8),  ph_corr_angle_vec.get(9),  ph_corr_angle_vec.get(10), ph_corr_angle_vec.get(11), 
        //        ph_corr_angle_vec.get(12), ph_corr_angle_vec.get(13), ph_corr_angle_vec.get(14), ph_corr_angle_vec.get(15));

        // Calculate the sin and cos of ph_corr_angle and store as a cfloat (cos in the real part, sin in the imaginary)
        auto ph_corr_vec = aie::sincos_complex(ph_corr_angle_vec);
        ph_corr_vec = aie::neg(ph_corr_vec);
        auto ph_corr_real = aie::real(ph_corr_vec);
        auto ph_corr_imag = aie::imag(ph_corr_vec);
        //printf("%d: ph_corr=[[%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], " \
        //       "[%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f]]\n", m_id, 
        //        ph_corr_real.get(0),  ph_corr_imag.get(0),   ph_corr_real.get(1),  ph_corr_imag.get(1),  
        //        ph_corr_real.get(2),  ph_corr_imag.get(2),   ph_corr_real.get(3),  ph_corr_imag.get(3),  
        //        ph_corr_real.get(4),  ph_corr_imag.get(4),   ph_corr_real.get(5),  ph_corr_imag.get(5),  
        //        ph_corr_real.get(6),  ph_corr_imag.get(6),   ph_corr_real.get(7),  ph_corr_imag.get(7),
        //        ph_corr_real.get(8),  ph_corr_imag.get(8),   ph_corr_real.get(9),  ph_corr_imag.get(9),  
        //        ph_corr_real.get(10), ph_corr_imag.get(10),  ph_corr_real.get(11), ph_corr_imag.get(11),  
        //        ph_corr_real.get(12), ph_corr_imag.get(12),  ph_corr_real.get(13), ph_corr_imag.get(13),  
        //        ph_corr_real.get(14), ph_corr_imag.get(14),  ph_corr_real.get(15), ph_corr_imag.get(15));
        
        
        //aie::maxdiff()

        //auto target_px_st = ph_corr_vec[0];
        //auto target_px_en = ph_corr_vec[15];
        


        //low_range_bins = aie::sub(low_range_bins, fft_seg);
        //auto high_range_bins = aie::add(low_range_bins, 1);

        //auto rc_vec[i] = rc_in_iter[range_bins[i]-(m_id*2048)];
        //printf("rc_vec[%d] = {%f, %f}\n", i, ((cfloat)rc_vec[i]).real, ((cfloat)rc_vec[i]).imag);


        //printf("m_range_grid.samp_start_size: %f\n", m_range_grid.samp_start_size);
        
        
        // Calculate the approximate index for differ_range_vec in the equally spaced range grid
        auto px_idx_acc = aie::mul(m_range_grid.inv_range_res, differ_range_vec);
        
        // Shift indices so they are not negative
        px_idx_acc = aie::add(px_idx_acc, 4096.0f);//-m_range_grid.seg_offset);
        //if (m_id == 0)
        //    px_idx_acc = aie::sub(px_idx_acc, 6144.0f);
        //else if (m_id == 1)
        //    px_idx_acc = aie::sub(px_idx_acc, 4096.0f);
        //else if (m_id == 2)
        //    px_idx_acc = aie::sub(px_idx_acc, 2048.0f);


        //auto px_idx_vec = aie::mul(m_range_grid.inv_range_res, differ_range_vec).to_vector<float>(0);
        //printf("m_range_grid.inv_range_res: %f\n", m_range_grid.inv_range_res);
        //print_float("px_idx_vec", px_idx_vec);
        //printf("%d: px_idx_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //       px_idx_vec.get(0),  px_idx_vec.get(1),  px_idx_vec.get(2),  px_idx_vec.get(3),
        //       px_idx_vec.get(4),  px_idx_vec.get(5),  px_idx_vec.get(6),  px_idx_vec.get(7),
        //       px_idx_vec.get(8),  px_idx_vec.get(9),  px_idx_vec.get(10), px_idx_vec.get(11),
        //       px_idx_vec.get(12), px_idx_vec.get(13), px_idx_vec.get(14), px_idx_vec.get(15));
        
        // Offset by samp_start_size (accommodates for fftshift)
        //px_idx_acc = aie::add(px_idx_acc, m_range_grid.samp_start_size);
        
        auto low_idx_float_vec = aie::sub(px_idx_acc, 0.5f).to_vector<float>(0);
        //auto high_idx_float_vec = aie::add(px_idx_acc, 0.5f).to_vector<float>(0);
        //printf("%d: low_idx_float_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //       low_idx_float_vec.get(0),  low_idx_float_vec.get(1),  low_idx_float_vec.get(2),  low_idx_float_vec.get(3),
        //       low_idx_float_vec.get(4),  low_idx_float_vec.get(5),  low_idx_float_vec.get(6),  low_idx_float_vec.get(7),
        //       low_idx_float_vec.get(8),  low_idx_float_vec.get(9),  low_idx_float_vec.get(10), low_idx_float_vec.get(11),
        //       low_idx_float_vec.get(12), low_idx_float_vec.get(13), low_idx_float_vec.get(14), low_idx_float_vec.get(15));

        // Round to nearest whole number
        auto low_idx_int_vec = aie::to_fixed<int32>(low_idx_float_vec);
        auto high_idx_int_vec = aie::add(low_idx_int_vec, 1);
        low_idx_float_vec = aie::to_float(low_idx_int_vec);
        print_int32("low_idx_int_vec", low_idx_int_vec);
        print_int32("high_idx_int_vec", high_idx_int_vec);
        //print_float("low_idx_float_vec", low_idx_float_vec);

        // Fractional part for interpolation
        auto px_delta_idx_vec = aie::sub(px_idx_acc, low_idx_float_vec).to_vector<float>(0);
        //printf("%d: px_delta_idx_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //       px_delta_idx_vec.get(0),  px_delta_idx_vec.get(1),  px_delta_idx_vec.get(2),  px_delta_idx_vec.get(3),
        //       px_delta_idx_vec.get(4),  px_delta_idx_vec.get(5),  px_delta_idx_vec.get(6),  px_delta_idx_vec.get(7),
        //       px_delta_idx_vec.get(8),  px_delta_idx_vec.get(9),  px_delta_idx_vec.get(10), px_delta_idx_vec.get(11),
        //       px_delta_idx_vec.get(12), px_delta_idx_vec.get(13), px_delta_idx_vec.get(14), px_delta_idx_vec.get(15));

        
        aie::vector<cfloat,16> interp_vec;
        //aie::vector<float,16> tmp2;
        //aie::vector<float,16> tmp3;
        printf("%d: rc_in_iter[0].real=%f\n", m_id, rc_in_iter[0].real);
        for(int i=0; i<16; i++) {
            //printf("%d: low_idx_int_vec[%d] = %d | high_idx_int_vec[%d] = %d\n", m_id, i, low_idx_int_vec.get(i), i, high_idx_int_vec.get(i));
            printf("%d: low_idx_int_vec[%d] = %d | rc_in_iter[low_idx_int_vec.get(i)].real = %f\n", m_id, i, low_idx_int_vec.get(i), rc_in_iter[low_idx_int_vec.get(i)].real);
            //printf("%d: rc_in_iter[low]=[%f, %f] | rc_in_iter[high]=[%f, %f]\n", m_id, 
            //        rc_in_iter[low_idx_int_vec.get(i)].real, rc_in_iter[low_idx_int_vec.get(i)].imag, 
            //        rc_in_iter[high_idx_int_vec.get(i)].real, rc_in_iter[high_idx_int_vec.get(i)].imag);
            auto rc_delta = aie::sub(rc_in_iter[high_idx_int_vec.get(i)], rc_in_iter[low_idx_int_vec.get(i)]);
            ////printf("%d: rc_delta=[%f, %f]\n", m_id, rc_delta.real, rc_delta.imag);
            auto px_rc_delta = aie::mul(rc_delta, (float)(px_delta_idx_vec.get(i)));
            interp_vec[i] = aie::add(px_rc_delta, rc_in_iter[low_idx_int_vec[i]]);
            ////tmp2[i] = rc_in_iter[low_idx_int_vec[i]].real + px_delta_idx_vec[i] * (rc_in_iter[high_idx_int_vec[i]].real - rc_in_iter[low_idx_int_vec[i]].real);
            ////tmp3[i] = rc_in_iter[low_idx_int_vec[i]].imag + px_delta_idx_vec[i] * (rc_in_iter[high_idx_int_vec[i]].imag - rc_in_iter[low_idx_int_vec[i]].imag);
        }
        //auto img_vec = aie::mul(interp_vec, ph_corr_vec).to_vector<cfloat>(0);
        auto img_real_vec = aie::real(interp_vec);
        auto img_imag_vec = aie::imag(interp_vec);
        printf("%d: img_vec=[[%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f]\n", m_id,
               //"[%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f]]\n", m_id, 
                img_real_vec.get(0),  img_imag_vec.get(0),   img_real_vec.get(1),  img_imag_vec.get(1),  
                img_real_vec.get(2),  img_imag_vec.get(2),   img_real_vec.get(3),  img_imag_vec.get(3),  
                img_real_vec.get(4),  img_imag_vec.get(4),   img_real_vec.get(5),  img_imag_vec.get(5),  
                img_real_vec.get(6),  img_imag_vec.get(6),   img_real_vec.get(7),  img_imag_vec.get(7));
                //img_real_vec.get(8),  img_imag_vec.get(8),   img_real_vec.get(9),  img_imag_vec.get(9),  
                //img_real_vec.get(10), img_imag_vec.get(10),  img_real_vec.get(11), img_imag_vec.get(11),  
                //img_real_vec.get(12), img_imag_vec.get(12),  img_real_vec.get(13), img_imag_vec.get(13),  
                //img_real_vec.get(14), img_imag_vec.get(14),  img_real_vec.get(15), img_imag_vec.get(15));
        //printf("%d: tmp2=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //       tmp2.get(0),  tmp2.get(1),  tmp2.get(2),  tmp2.get(3),
        //       tmp2.get(4),  tmp2.get(5),  tmp2.get(6),  tmp2.get(7),
        //       tmp2.get(8),  tmp2.get(9),  tmp2.get(10), tmp2.get(11),
        //       tmp2.get(12), tmp2.get(13), tmp2.get(14), tmp2.get(15));
        //printf("%d: tmp3=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //       tmp3.get(0),  tmp3.get(1),  tmp3.get(2),  tmp3.get(3),
        //       tmp3.get(4),  tmp3.get(5),  tmp3.get(6),  tmp3.get(7),
        //       tmp3.get(8),  tmp3.get(9),  tmp3.get(10), tmp3.get(11),
        //       tmp3.get(12), tmp3.get(13), tmp3.get(14), tmp3.get(15));


        //for(int i=0; i<2048/16; i++) {
        //    auto rc_vec = *rc_in_iter++; // Answer; interp or nearest neighbors

        //    for(int j=0; j<16; j++) {
        //        float range_res = (j + i*16)*m_range_grid.range_res; // Actual real x values
        //        auto target_px = differ_range_vec[j]; //ranges of target X


        //        if (target_px > rc_vec[0] && target_px < rc_vec[15]);
        //    }
        //}

        auto prev_img_vec = aie::load_v<16>(m_img + xy_idx*16);
        auto updated_img_vec = aie::add(ph_corr_vec, prev_img_vec);
        updated_img_vec.store(m_img + xy_idx*16);

        //printf("%d: ph_corr=[[%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f]]\n", m_id,
        //        ph_corr_real.get(0),  ph_corr_imag.get(0),   ph_corr_real.get(1),  ph_corr_imag.get(1),  
        //        ph_corr_real.get(2),  ph_corr_imag.get(2),   ph_corr_real.get(3),  ph_corr_imag.get(3));
        
        //// Determine which pixels fall within the swath range
        ////auto right_mask = aie::ge(differ_range, min_range_bins);
        ////auto left_mask = aie::le(differ_range, max_range_bins);
        ////auto differ_range_mask = right_mask & left_mask;
        //
        //auto px_diff = aie::add(m_range_grid.range_width/2, x_pxls);
        //printf("%d: px_diff=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
        //        px_diff.get(0),  px_diff.get(1),  px_diff.get(2),  px_diff.get(3), 
        //        px_diff.get(4),  px_diff.get(5),  px_diff.get(6),  px_diff.get(7),
        //        px_diff.get(8),  px_diff.get(9),  px_diff.get(10), px_diff.get(11), 
        //        px_diff.get(12), px_diff.get(13), px_diff.get(14), px_diff.get(15));

        ////aie::saturation_mode current_sat=aie::get_saturation();
        ////printf("CURRENT_SATURATION = %d\n", current_sat);

        //auto range_bins_float = aie::mul(px_diff, m_range_grid.inv_range_res).to_vector<float>(0);
        //printf("%d: range_bins_float=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
        //        range_bins_float.get(0),  range_bins_float.get(1),  range_bins_float.get(2),  range_bins_float.get(3), 
        //        range_bins_float.get(4),  range_bins_float.get(5),  range_bins_float.get(6),  range_bins_float.get(7),
        //        range_bins_float.get(8),  range_bins_float.get(9),  range_bins_float.get(10), range_bins_float.get(11), 
        //        range_bins_float.get(12), range_bins_float.get(13), range_bins_float.get(14), range_bins_float.get(15));

        //// TODO: This is actually rounding correctly...but I need it to floor round...
        //aie::vector<int32,16> range_bins = aie::to_fixed<int32>(range_bins_float);
        //printf("%d: range_bins=%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", m_id, 
        //        range_bins.get(0),  range_bins.get(1),  range_bins.get(2),  range_bins.get(3), 
        //        range_bins.get(4),  range_bins.get(5),  range_bins.get(6),  range_bins.get(7),
        //        range_bins.get(8),  range_bins.get(9),  range_bins.get(10), range_bins.get(11), 
        //        range_bins.get(12), range_bins.get(13), range_bins.get(14), range_bins.get(15));

        //low_range_bins = aie::sub(low_range_bins, fft_seg);
        //auto high_range_bins = aie::add(low_range_bins, 1);

        //auto rc_vec[i] = rc_in_iter[range_bins[i]-(m_id*2048)];
        //printf("rc_vec[%d] = {%f, %f}\n", i, ((cfloat)rc_vec[i]).real, ((cfloat)rc_vec[i]).imag);

        //y0+(x-x0)*((y1 - y0)/dr);
        // Interpolation...HOW TO DO THIS
        //y0+(x-x0)*((y1 - y0)/(x1-x0));
        
        //if (differ_range_mask.full()) {
        //    auto img_vec = aie::mul(differ_range, ph_corr).to_vector<cfloat>(0);
        //    aie::store_v(img+(y_idx*x_px_grid_len)+(x_idx*16), img_vec);
        //}

        //// Increment range bins to next segment
        //range_bins = aie::add(dr, range_bins);
    }

    m_iter++;
    
    // Check if image output should sent
    if (m_iter == ACCUM_PULSES) {
        auto img_iter = aie::begin_vector<16>(m_img);
        for(unsigned i=0; i<2; i++) {
            *img_out_iter++ = *img_iter++;
        }
    }

    img_out.release();
}
