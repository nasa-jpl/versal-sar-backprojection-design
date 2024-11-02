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
{}

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
}

void Backprojection::init_pixel_grid(float x_st, float x_en, int x_len, float y_st, float y_en, int y_len)
{
    m_px_grid.x_st = x_st;
    m_px_grid.x_en = x_en;
    m_px_grid.x_len = x_len;
    m_px_grid.x_res = (x_en - x_st) / (x_len - 1);
    m_px_grid.y_st = y_st;
    m_px_grid.y_en = y_en;
    m_px_grid.y_len = y_len;
    m_px_grid.y_res = (y_en - y_st) / (y_len - 1);
}

void Backprojection::init_pixel_segment()
{
    m_seg.range_st = ((m_range_grid.range_width/BP_SOLVERS)*m_id) - m_range_grid.range_width/2;
    m_seg.range_en = ((m_range_grid.range_width/BP_SOLVERS)*(m_id+1)) - m_range_grid.range_width/2;

    m_seg.st_x_px_bound = false;
    m_seg.en_x_px_bound = false;
    m_seg.all_x_px_bound = false;
    if (m_px_grid.x_st >= m_seg.range_st && m_px_grid.x_st < m_seg.range_en) {
        m_seg.st_x_px_bound = true;
    }
    if (m_px_grid.x_en >= m_seg.range_st && m_px_grid.x_en < m_seg.range_en) {
        m_seg.en_x_px_bound = true;
    }
    if (m_px_grid.x_st <= m_seg.range_st && m_px_grid.x_en > m_seg.range_en) {
        m_seg.all_x_px_bound = true;
    }
}

void Backprojection::backprojection_kern(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
                                         input_buffer<TT_DATA, extents<2048>>& __restrict rc_in,
                                         input_buffer<TT_DATA, extents<2048>>& __restrict xy_px_in,
                                         output_async_buffer<TT_DATA, extents<2048>>& __restrict img_out) {
    // Lock img_out buffer
    img_out.acquire();
    printf("%d: AUSTIN\n", m_id);

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

    // TODO: START USING FUNCTIONS; USE BETTER NAMING CONVENTIONS

    // Image memory
    //alignas(aie::vector_decl_align) TT_DATA img[m_px_grid.x_len*m_px_grid.y_len];
    

    // Initialize vector of pixels to iterate though
    //alignas(aie::vector_decl_align) float start_x_pxls[16];
    ////m_seg.x_pxls = aie::zeros<float,16>();
    ////m_seg.y_pxls = aie::zeros<float,16>();
    //aie::vector<float,16> x_pxls = aie::zeros<float,16>();
    //aie::vector<float,16> y_pxls = aie::zeros<float,16>();
    //m_seg.valid_elems = 0;
    //if (m_seg.st_x_px_bound) {
    //    float x_pxl;
    //    for(int i = 0; i < 16; i++) chess_prepare_for_pipelining {
    //        x_pxl = m_px_grid.x_st+(m_px_grid.x_res*i);
    //        if (x_pxl < m_seg.range_en) {
    //            x_pxls[i] = x_pxl;
    //            m_seg.valid_elems = i+1;
    //        }
    //    }
    //    for(int i = 0; i < 16; i++) chess_prepare_for_pipelining {
    //        y_pxls[i] = m_px_grid.y_st+(m_px_grid.y_res*i);
    //    }
    //}
    //else if (m_seg.en_x_px_bound || m_seg.all_x_px_bound) {
    //    // Calculate number of steps for new start position based on x_res (rounds up)
    //    int steps = (int)((m_seg.range_st-m_px_grid.x_st)/m_px_grid.x_res) + 1;
    //    float new_x_px_st = m_px_grid.x_st + steps*m_px_grid.x_res;
    //    float x_pxl;
    //    for(int i = 0; i < 16; i++) chess_prepare_for_pipelining {
    //        x_pxl = new_x_px_st+(m_px_grid.x_res*i);
    //        if (x_pxl < m_seg.range_en) {
    //            x_pxls[i] = x_pxl;
    //            m_seg.valid_elems = i+1;
    //        }
    //    }
    //    for(int i = 0; i < 16; i++) chess_prepare_for_pipelining {
    //        y_pxls[i] = m_px_grid.y_st+(m_px_grid.y_res*i);
    //    }
    //}
    //aie::store_v(start_x_pxls, x_pxls);
    //
    //// Increment x_pxls and y_pxls by this amount every iteration
    //float x_px_inc = m_px_grid.x_res*16;
    //float y_px_inc = m_px_grid.y_res*16;
    
    //if (m_seg.st_x_px_bound || m_seg.en_x_px_bound || m_seg.all_x_px_bound) {
    //    m_px_grid.display(m_id);
    //    m_seg.display(m_id);
    //    m_range_grid.display(m_id);
    //    printf("%d: x_pxls=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
    //            x_pxls.get(0),  x_pxls.get(1),  x_pxls.get(2),  x_pxls.get(3), 
    //            x_pxls.get(4),  x_pxls.get(5),  x_pxls.get(6),  x_pxls.get(7),
    //            x_pxls.get(8),  x_pxls.get(9),  x_pxls.get(10), x_pxls.get(11), 
    //            x_pxls.get(12), x_pxls.get(13), x_pxls.get(14), x_pxls.get(15));
    //    printf("%d: y_pxls=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
    //            y_pxls.get(0),  y_pxls.get(1),  y_pxls.get(2),  y_pxls.get(3), 
    //            y_pxls.get(4),  y_pxls.get(5),  y_pxls.get(6),  y_pxls.get(7),
    //            y_pxls.get(8),  y_pxls.get(9),  y_pxls.get(10), y_pxls.get(11), 
    //            y_pxls.get(12), y_pxls.get(13), y_pxls.get(14), y_pxls.get(15));
    //}
    //float dr = max_px_width/TP_POINT_SIZE;
    
    // Half of the sample size (used for shifting)
    //float half_fft_sz = TP_POINT_SIZE/2;

    // Min and max range bins
    //float min_range_bins = -(TP_POINT_SIZE/2)*dr;
    //float max_range_bins = ((TP_POINT_SIZE/2) - 1)*dr;
    //
    //// Initialize range bins
    //aie::vector<float,16> range_bins = aie::zeros<float,16>();
    //for(int i = 0; i < 16; i++) chess_prepare_for_pipelining {
    //    range_bins[i] = i;
    //}
    //
    //// Shift range by half_fft_sz to align with fftshifting in graph.h
    //range_bins = aie::sub(range_bins, half_fft_sz);

    //// Normalize and scale range bins to range distance
    //range_bins = aie::mul(dr, range_bins).to_vector<float>(0);

    // Precalculate z_diff_sq
    auto z_diff_sq = z_ant*z_ant;
    //auto a = aie::cos(1 << 1);
    //auto b = aie::cos(1 << 2);
    //auto j = aie::cos(1000000000);
    //auto jj = aie::cos(0x3FFFF000);
    //auto jjj = aie::cos(0x3FFFFC00);
    //auto jjjj = aie::cos(0x3FFFFFF0);
    //auto jjjjj = aie::cos(0x3FFFFFFF);
    //auto k = aie::cos(2000000000);
    //auto l = aie::cos(2147483647);
    //auto m = aie::cos(2147483648);
    //printf("%d: %d %d %d %d %d %d %d %d, %d, %d, %d, %d, %d, %d, %d, %d, %d %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d\n", m_id,a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,s,t,u,v,w,x,y,z,aa,bb,cc,dd,ee,ff);
    //for(int xy_idx=0; xy_idx < 2048/16; xy_idx++) chess_prepare_for_pipelining {
    for(int xy_idx=0; xy_idx < 1; xy_idx++) chess_prepare_for_pipelining {
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
        printf("%d: differ_range_vec=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
                differ_range_vec.get(0),  differ_range_vec.get(1),  differ_range_vec.get(2),  differ_range_vec.get(3), 
                differ_range_vec.get(4),  differ_range_vec.get(5),  differ_range_vec.get(6),  differ_range_vec.get(7),
                differ_range_vec.get(8),  differ_range_vec.get(9),  differ_range_vec.get(10), differ_range_vec.get(11), 
                differ_range_vec.get(12), differ_range_vec.get(13), differ_range_vec.get(14), differ_range_vec.get(15));

        // Calculate phase correction for image
        auto ph_corr_angle_vec = aie::mul(m_radar_params.ph_corr_coef, differ_range_vec).to_vector<float>(0);
        //aie::print(ph_corr_angle_vec, true, "ph_corr_angle_vec (before scale)=");
        printf("%d: ph_corr_angle_vec (before scale)=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
                ph_corr_angle_vec.get(0),  ph_corr_angle_vec.get(1),  ph_corr_angle_vec.get(2),  ph_corr_angle_vec.get(3), 
                ph_corr_angle_vec.get(4),  ph_corr_angle_vec.get(5),  ph_corr_angle_vec.get(6),  ph_corr_angle_vec.get(7),
                ph_corr_angle_vec.get(8),  ph_corr_angle_vec.get(9),  ph_corr_angle_vec.get(10), ph_corr_angle_vec.get(11), 
                ph_corr_angle_vec.get(12), ph_corr_angle_vec.get(13), ph_corr_angle_vec.get(14), ph_corr_angle_vec.get(15));

        // Floor round to neg infinity by casting to int32, then back to float for later operations
        auto num_pi_wrapped_acc = aie::mul(INV_TWO_PI, ph_corr_angle_vec);

        //auto neg_mask = aie::lt(num_pi_wrapped_vec, 0.0f);
        //aie::print(neg_mask, true, "neg_mask=");
        //printf("%d: neg_mask=[%d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d, %d]\n", m_id, 
        //        neg_mask.test(0),  neg_mask.test(1),  neg_mask.test(2),  neg_mask.test(3), 
        //        neg_mask.test(4),  neg_mask.test(5),  neg_mask.test(6),  neg_mask.test(7),
        //        neg_mask.test(8),  neg_mask.test(9),  neg_mask.test(10), neg_mask.test(11), 
        //        neg_mask.test(12), neg_mask.test(13), neg_mask.test(14), neg_mask.test(15));
        auto num_pi_wrapped_vec = aie::sub(num_pi_wrapped_acc, 0.5f).to_vector<float>(0);
        //auto num_pi_wrapped_floor_vec = aie::sub(num_pi_wrapped_vec, 0.5f);
        //auto num_pi_wrapped_ceil_vec = aie::add(num_pi_wrapped_vec, 0.5f);
        //num_pi_wrapped_vec = aie::select(num_pi_wrapped_floor_vec, num_pi_wrapped_ceil_vec, neg_mask);
        auto num_pi_wrapped_int_vec = aie::to_fixed<int32>(num_pi_wrapped_vec);
        num_pi_wrapped_vec = aie::to_float(num_pi_wrapped_int_vec);
        //aie::print(num_pi_wrapped_vec, true, "num_pi_wrapped_vec=");
        printf("%d: num_pi_wrapped_vec=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
                num_pi_wrapped_vec.get(0),  num_pi_wrapped_vec.get(1),  num_pi_wrapped_vec.get(2),  num_pi_wrapped_vec.get(3), 
                num_pi_wrapped_vec.get(4),  num_pi_wrapped_vec.get(5),  num_pi_wrapped_vec.get(6),  num_pi_wrapped_vec.get(7),
                num_pi_wrapped_vec.get(8),  num_pi_wrapped_vec.get(9),  num_pi_wrapped_vec.get(10), num_pi_wrapped_vec.get(11), 
                num_pi_wrapped_vec.get(12), num_pi_wrapped_vec.get(13), num_pi_wrapped_vec.get(14), num_pi_wrapped_vec.get(15));

        auto scale_down_angle_acc = aie::negmul(TWO_PI, num_pi_wrapped_vec);
        scale_down_angle_acc = aie::sub(scale_down_angle_acc, PI);
        aie::print(scale_down_angle_acc, true, "scale_down_angle_acc=");
        //printf("%d: scale_down_angle=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
        //        scale_down_angle_acc[0],  scale_down_angle[1],  scale_down_angle[2],  scale_down_angle[3], 
        //        scale_down_angle_acc[4],  scale_down_angle[5],  scale_down_angle[6],  scale_down_angle[7],
        //        scale_down_angle_acc[8],  scale_down_angle[9],  scale_down_angle[10], scale_down_angle[11], 
        //        scale_down_angle_acc[12], scale_down_angle[13], scale_down_angle[14], scale_down_angle[15]);
        ph_corr_angle_vec = aie::add(scale_down_angle_acc, ph_corr_angle_vec).to_vector<float>(0);

        //aie::print(ph_corr_angle_vec, true, "ph_corr_angle_vec (after scale)=");
        printf("%d: ph_corr_angle (after scale)=[%f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f, %f]\n", m_id, 
                ph_corr_angle_vec.get(0),  ph_corr_angle_vec.get(1),  ph_corr_angle_vec.get(2),  ph_corr_angle_vec.get(3), 
                ph_corr_angle_vec.get(4),  ph_corr_angle_vec.get(5),  ph_corr_angle_vec.get(6),  ph_corr_angle_vec.get(7),
                ph_corr_angle_vec.get(8),  ph_corr_angle_vec.get(9),  ph_corr_angle_vec.get(10), ph_corr_angle_vec.get(11), 
                ph_corr_angle_vec.get(12), ph_corr_angle_vec.get(13), ph_corr_angle_vec.get(14), ph_corr_angle_vec.get(15));
        auto ph_corr_vec = aie::sincos_complex(ph_corr_angle_vec);
        ph_corr_vec = aie::neg(ph_corr_vec);
        auto ph_corr_real = aie::real(ph_corr_vec);
        auto ph_corr_imag = aie::imag(ph_corr_vec);
        printf("%d: ph_corr=[[%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], " \
               "[%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f], [%f, %f]]\n", m_id, 
                ph_corr_real.get(0),  ph_corr_imag.get(0),   ph_corr_real.get(1),  ph_corr_imag.get(1),  
                ph_corr_real.get(2),  ph_corr_imag.get(2),   ph_corr_real.get(3),  ph_corr_imag.get(3),  
                ph_corr_real.get(4),  ph_corr_imag.get(4),   ph_corr_real.get(5),  ph_corr_imag.get(5),  
                ph_corr_real.get(6),  ph_corr_imag.get(6),   ph_corr_real.get(7),  ph_corr_imag.get(7),  
                ph_corr_real.get(8),  ph_corr_imag.get(8),   ph_corr_real.get(9),  ph_corr_imag.get(9),  
                ph_corr_real.get(10), ph_corr_imag.get(10),  ph_corr_real.get(11), ph_corr_imag.get(11),  
                ph_corr_real.get(12), ph_corr_imag.get(12),  ph_corr_real.get(13), ph_corr_imag.get(13),  
                ph_corr_real.get(14), ph_corr_imag.get(14),  ph_corr_real.get(15), ph_corr_imag.get(15));

        *img_out_iter++ = ph_corr_vec;

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
}

    //aie::vector<TT_DATA,16> rc_vec;
    ////printf("ID: %d\n", m_id);
    //if (m_seg.st_x_px_bound || m_seg.en_x_px_bound || m_seg.all_x_px_bound) {
    //    //printf("ID: %d\n", m_id);
    //    for(int y_idx=0; y_idx < m_px_grid.y_len/16; y_idx++) chess_prepare_for_pipelining {
    //        for(int x_idx=0; x_idx < m_px_grid.x_len/16; x_idx++) chess_prepare_for_pipelining {

    //            // Calculate differential range for pixel segments
    //            auto x = aie::sub(x_ant, x_pxls);
    //            auto x_diff_sq = aie::mul_square(x).to_vector<float>(0);

    //            auto y = aie::sub(y_ant, y_pxls);
    //            auto y_diff_sq = aie::mul_square(y).to_vector<float>(0);

    //            auto xy_add = aie::add(x_diff_sq, y_diff_sq);
    //            auto xyz_add = aie::add(xy_add, z_diff_sq);

    //            auto xyz_sqrt = aie::sqrt(xyz_add);

    //            auto differ_range = aie::sub(xyz_sqrt, r0);

    //            // Calculate phase correction for image
    //            auto ph_corr_angle = aie::mul(m_radar_params.ph_corr_coef, differ_range).to_vector<float>(0);
    //            auto ph_corr = aie::sincos_complex(ph_corr_angle);
    //            
    //            // Determine which pixels fall within the swath range
    //            //auto right_mask = aie::ge(differ_range, min_range_bins);
    //            //auto left_mask = aie::le(differ_range, max_range_bins);
    //            //auto differ_range_mask = right_mask & left_mask;
    //            
    //            //TODO: WRONG
    //            auto px_diff = aie::add(m_range_grid.range_width/2, x_pxls);
    //            //printf("%d: px_diff=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
    //            //        px_diff.get(0),  px_diff.get(1),  px_diff.get(2),  px_diff.get(3), 
    //            //        px_diff.get(4),  px_diff.get(5),  px_diff.get(6),  px_diff.get(7),
    //            //        px_diff.get(8),  px_diff.get(9),  px_diff.get(10), px_diff.get(11), 
    //            //        px_diff.get(12), px_diff.get(13), px_diff.get(14), px_diff.get(15));

    //            //aie::saturation_mode current_sat=aie::get_saturation();
    //            //printf("CURRENT_SATURATION = %d\n", current_sat);

    //            auto range_bins_float = aie::mul(px_diff, m_range_grid.inv_range_res).to_vector<float>(0);
    //            printf("%d: range_bins_float=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
    //                    range_bins_float.get(0),  range_bins_float.get(1),  range_bins_float.get(2),  range_bins_float.get(3), 
    //                    range_bins_float.get(4),  range_bins_float.get(5),  range_bins_float.get(6),  range_bins_float.get(7),
    //                    range_bins_float.get(8),  range_bins_float.get(9),  range_bins_float.get(10), range_bins_float.get(11), 
    //                    range_bins_float.get(12), range_bins_float.get(13), range_bins_float.get(14), range_bins_float.get(15));

    //            // TODO: This is actually rounding correctly...but I need it to floor round...
    //            aie::vector<int32,16> range_bins = aie::to_fixed<int32>(range_bins_float);
    //            printf("%d: range_bins=%d %d %d %d %d %d %d %d %d %d %d %d %d %d %d %d\n", m_id, 
    //                    range_bins.get(0),  range_bins.get(1),  range_bins.get(2),  range_bins.get(3), 
    //                    range_bins.get(4),  range_bins.get(5),  range_bins.get(6),  range_bins.get(7),
    //                    range_bins.get(8),  range_bins.get(9),  range_bins.get(10), range_bins.get(11), 
    //                    range_bins.get(12), range_bins.get(13), range_bins.get(14), range_bins.get(15));

    //            //low_range_bins = aie::sub(low_range_bins, fft_seg);
    //            //auto high_range_bins = aie::add(low_range_bins, 1);

    //            for (int i=0; i<m_seg.valid_elems; i++) {
    //                rc_vec[i] = rc_in_iter[range_bins[i]-(m_id*2048)];
    //                printf("rc_vec[%d] = {%f, %f}\n", i, ((cfloat)rc_vec[i]).real, ((cfloat)rc_vec[i]).imag);
    //            }

    //            //y0+(x-x0)*((y1 - y0)/dr);
    //            // Interpolation...HOW TO DO THIS
    //            //y0+(x-x0)*((y1 - y0)/(x1-x0));
    //            
    //            //if (differ_range_mask.full()) {
    //            //    auto img_vec = aie::mul(differ_range, ph_corr).to_vector<cfloat>(0);
    //            //    aie::store_v(img+(y_idx*x_px_grid_len)+(x_idx*16), img_vec);
    //            //}

    //            //// Increment range bins to next segment
    //            //range_bins = aie::add(dr, range_bins);
    //        }
    //    }
    //}
//}

    //img_out.release();
    
    //aie::vector<TT_DATA,16> low_ph_data;
    //aie::vector<TT_DATA,16> high_ph_data;
    //if (m_seg.st_x_px_bound || m_seg.en_x_px_bound || m_seg.all_x_px_bound) {
    //    //printf("ID: %d\n", m_id);
    //    for(int y_idx=0; y_idx < m_px_grid.y_len/16; y_idx++) chess_prepare_for_pipelining {
    //        for(int x_idx=0; x_idx < m_px_grid.x_len/16; x_idx++) chess_prepare_for_pipelining {

    //            // Calculate differential range for pixel segments
    //            auto x = aie::sub(x_ant, x_pxls);
    //            auto x_diff_sq = aie::mul_square(x).to_vector<float>(0);

    //            auto y = aie::sub(y_ant, y_pxls);
    //            auto y_diff_sq = aie::mul_square(y).to_vector<float>(0);

    //            auto xy_add = aie::add(x_diff_sq, y_diff_sq);
    //            auto xyz_add = aie::add(xy_add, z_diff_sq);

    //            auto xyz_sqrt = aie::sqrt(xyz_add);

    //            auto differ_range = aie::sub(xyz_sqrt, r0);

    //            // Calculate phase correction for image
    //            auto ph_corr_angle = aie::mul(m_radar_params.ph_corr_coef, differ_range).to_vector<float>(0);
    //            auto ph_corr = aie::sincos_complex(ph_corr_angle);
    //            
    //            // Determine which pixels fall within the swath range
    //            //auto right_mask = aie::ge(differ_range, min_range_bins);
    //            //auto left_mask = aie::le(differ_range, max_range_bins);
    //            //auto differ_range_mask = right_mask & left_mask;
    //            
    //            //TODO: WRONG
    //            auto px_diff = aie::sub(m_range_grid.range_width, x_pxls);
    //            printf("%d: px_diff=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
    //                    px_diff.get(0),  px_diff.get(1),  px_diff.get(2),  px_diff.get(3), 
    //                    px_diff.get(4),  px_diff.get(5),  px_diff.get(6),  px_diff.get(7),
    //                    px_diff.get(8),  px_diff.get(9),  px_diff.get(10), px_diff.get(11), 
    //                    px_diff.get(12), px_diff.get(13), px_diff.get(14), px_diff.get(15));

    //            //aie::saturation_mode current_sat=aie::get_saturation();
    //            //printf("CURRENT_SATURATION = %d\n", current_sat);

    //            auto low_range_bins_float = aie::mul(px_diff, m_range_grid.inv_range_res).to_vector<float>(0);
    //            printf("%d: low_range_bins_float=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
    //                    low_range_bins_float.get(0),  low_range_bins_float.get(1),  low_range_bins_float.get(2),  low_range_bins_float.get(3), 
    //                    low_range_bins_float.get(4),  low_range_bins_float.get(5),  low_range_bins_float.get(6),  low_range_bins_float.get(7),
    //                    low_range_bins_float.get(8),  low_range_bins_float.get(9),  low_range_bins_float.get(10), low_range_bins_float.get(11), 
    //                    low_range_bins_float.get(12), low_range_bins_float.get(13), low_range_bins_float.get(14), low_range_bins_float.get(15));

    //            // TODO: This is actually rounding correctly...but I need it to floor round...
    //            aie::vector<int32,16> low_range_bins = aie::to_fixed<int32>(low_range_bins_float);
    //            printf("%d: low_range_bins=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
    //                    low_range_bins.get(0),  low_range_bins.get(1),  low_range_bins.get(2),  low_range_bins.get(3), 
    //                    low_range_bins.get(4),  low_range_bins.get(5),  low_range_bins.get(6),  low_range_bins.get(7),
    //                    low_range_bins.get(8),  low_range_bins.get(9),  low_range_bins.get(10), low_range_bins.get(11), 
    //                    low_range_bins.get(12), low_range_bins.get(13), low_range_bins.get(14), low_range_bins.get(15));

    //            //low_range_bins = aie::sub(low_range_bins, fft_seg);
    //            //auto high_range_bins = aie::add(low_range_bins, 1);

    //            
    //            //for (int i=0; i<16; i++) {
    //            //    low_ph_data[i] = ph_data_iter[low_range_bins[i]];
    //            //    high_ph_data[i] = ph_data_iter[low_range_bins[i]+1];
    //            //    printf("ph_corr[%d] = {%f, %f}\n", i, ((cfloat)ph_corr[i]).real, ((cfloat)ph_corr[i]).imag);
    //            //    printf("low_ph_data[%d] = {%f, %f}\n", i, ((cfloat)low_ph_data[i]).real, ((cfloat)low_ph_data[i]).imag);
    //            //    printf("high_ph_data[%d] = {%f, %f}\n", i, ((cfloat)high_ph_data[i]).real, ((cfloat)high_ph_data[i]).imag);
    //            //}

    //            //y0+(x-x0)*((y1 - y0)/dr);
    //            // Interpolation...HOW TO DO THIS
    //            //y0+(x-x0)*((y1 - y0)/(x1-x0));
    //            
    //            //if (differ_range_mask.full()) {
    //            //    auto img_vec = aie::mul(differ_range, ph_corr).to_vector<cfloat>(0);
    //            //    aie::store_v(img+(y_idx*x_px_grid_len)+(x_idx*16), img_vec);
    //            //}

    //            //// Increment range bins to next segment
    //            //range_bins = aie::add(dr, range_bins);
    //            
    //            // Increment X pixel segment
    //            x_pxls = aie::add(x_pxls, x_px_inc);
    //        }

    //        // Reset X pixel segment
    //        x_pxls = aie::load_v<16>(start_x_pxls);

    //        // Increment Y pixel segment
    //        y_pxls = aie::add(y_pxls, y_px_inc);
    //    }
    //}

    ////img_out.release();
//}

//void phase_corr_kern(input_buffer<float, extents<2048>> &restrict x_ant_pos,
//                     input_buffer<float, extents<2048>> &restrict y_ant_pos,
//                     input_buffer<float, extents<2048>> &restrict z_ant_pos,
//                     input_buffer<float, extents<2048>> &restrict ref_range,
//                     input_buffer<TT_DATA, extents<2048>> &restrict ph_data,
//                     output_buffer<float, extents<2048>> &restrict img_out) {
//    
//    // Antenna X, Y, and Z position
//    auto x_ant_pos_iter = aie::begin(x_ant_pos);
//    auto y_ant_pos_iter = aie::begin(y_ant_pos);
//    auto z_ant_pos_iter = aie::begin(z_ant_pos);
//
//    // Reference range to scene center
//    auto ref_range_iter = aie::begin(ref_range);
//
//    // Phase history data in time domain (compressed range lines)
//    auto ph_data_circ_iter = aie::cbegin_vector_circular<16, 2048>(ph_data.data());
//
//    // Image output
//    auto img_out_iter = aie::begin_vector<16>(img_out);
//
//    // Pixel X start and end position on grid
//    float x_px_st = -5.0;
//    float x_px_en = 5.0;
//
//    // Pixel Y start and end position on grid
//    float y_px_st = -5.0;
//    float y_px_en = 5.0;
//    
//    // Pixel X & Y grid length (128x128)
//    int x_px_grid_len = 128;
//    int y_px_grid_len = 128;
//    
//    // Pixel X & Y step sizes
//    float dx = (x_px_en - x_px_st) / (x_px_grid_len - 1);
//    float dy = (y_px_en - y_px_st) / (y_px_grid_len - 1);
//
//    // Initialize vector of pixels to iterate though (maybe there is 
//    // a more efficient way to do this with circular buffers?)
//    aie::vector<float,32> x_pxls = aie::zeros<float,32>();
//    aie::vector<float,32> y_pxls = aie::zeros<float,32>();
//    for(int i = 0; i < 32; i++) chess_prepare_for_pipelining {
//        x_pxls[i] = x_px_st+(dx*i);
//        y_pxls[i] = y_px_st+(dy*i);
//    }
//    
//    // Increment x_pxls and y_pxls by this amount every iteration
//    aie::vector<float,32> x_px_inc = aie::broadcast<float,32>(dx*32);
//    aie::vector<float,32> y_px_inc = aie::broadcast<float,32>(dy*32);
//
//    // Extract antenna position and range to scene center for each pulse
//    float x_ant = *x_ant_pos_iter++;
//    float y_ant = *y_ant_pos_iter++;
//    float z_ant = *z_ant_pos_iter++;
//    float r0 = *ref_range_iter++;
//    
//    //Constants
//    const float pi = 3.141593565;
//    const float minF = 100000.0;
//    const float c = 299792458.0;
//    const int dF = 10000;
//    const int sample_size = 8192;
//
//    // Calculate coefficient in phase correction
//    auto ph_corr_coef = (4*pi*minF)/c;
//
//    // Calculate max width scene size of image (in range direction)
//    auto max_px_width = c/(2*dF);
//    
//    // Precalculate z_diff_sq
//    auto z_diff_sq = z_ant*z_ant;
//    
//    int x_px_idx, y_px_idx;
//    for(unsigned p = 0; p < 2048; p++) chess_prepare_for_pipelining {
//        for(y_px_idx = 0; y_px_idx < y_px_grid_len/32; y_px_idx++) chess_prepare_for_pipelining {
//            for(x_px_idx = 0; x_px_idx < x_px_grid_len/32; x_px_idx++) chess_prepare_for_pipelining {
//
//                // Calculate differential range for pixel segments
//                auto x = aie::sub(x_ant, x_pxls);
//                auto x_diff_sq = aie::mul_square(x).to_vector<float>(0);
//
//                auto y = aie::sub(y_ant, y_pxls);
//                auto y_diff_sq = aie::mul_square(y).to_vector<float>(0);
//
//                auto xy_add = aie::add(x_diff_sq, y_diff_sq);
//                auto xyz_add = aie::add(xy_add, z_diff_sq);
//
//                auto xyz_sqrt = aie::sqrt(xyz_add);
//
//                auto dR = aie::sub(xyz_sqrt, r0);
//
//                // Calculate phase correction for image
//                auto ph_corr_angle = aie::mul(ph_corr_coef, dR).to_vector<float>(0);
//                auto [ph_corr_imag, ph_corr_real] = aie::sincos(ph_corr_angle);
//                
//                // Doesn't work. Can't store a vector<cflaot, 32> :(
//                //auto ph_corr = aie::sincos_complex(ph_corr_angle);
//                
//                //aie::sub(dR-)
//
//                //y0+(x-x0)*((y1 - y0)/(x1-x0));
//                
//                // Increment X pixel segment
//                x_pxls = aie::add(x_pxls, x_px_inc);
//            }
//
//            // Reset X pixel segment
//            auto x_px_reset = aie::mul((float)x_px_idx, x_px_inc).to_vector<float>(0);
//            x_pxls = aie::sub(x_pxls, x_px_reset);
//
//            // Increment Y pixel segment
//            y_pxls = aie::add(y_pxls, y_px_inc);
//        }
//
//        // Reset X and Y pixel segments
//        auto x_px_reset = aie::mul((float)x_px_idx, x_px_inc).to_vector<float>(0);
//        auto y_px_reset = aie::mul((float)y_px_idx, y_px_inc).to_vector<float>(0);
//
//        x_pxls = aie::sub(x_pxls, x_px_reset);
//        y_pxls = aie::sub(y_pxls, y_px_reset);
//        
//        // Increment to next range line data
//        x_ant = *x_ant_pos_iter++;
//        y_ant = *y_ant_pos_iter++;
//        z_ant = *z_ant_pos_iter++;
//        r0 = *ref_range_iter++;
//        z_diff_sq = z_ant*z_ant;
//
//        //auto data = aie::conj(*in_iter++);
//        //printf("i: %d\n", i);
//        //aie::print(data, true, "Data:\n");
//        //*out_iter++ = data;
//        //
//        //y0+(x-x0)*((y1 - y0)/(x1-x0));
//        //*img_out_iter++ = aie::abs();
//    }
//}
