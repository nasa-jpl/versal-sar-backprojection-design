// By: Austin Owens
// Date: 10/16/2024
// Desc: Performs phase correction compressed range lines for target

#include "aie_api/aie_adf.hpp"
#include "aie_api/utils.hpp"
#include "custom_kernels.h"

using namespace adf;

void slowtime_splicer_kern(input_buffer<float, extents<1>>& __restrict x_ant_pos,
                           input_buffer<float, extents<1>>& __restrict y_ant_pos,
                           input_buffer<float, extents<1>>& __restrict z_ant_pos,
                           input_buffer<float, extents<1>>& __restrict ref_range,
                           output_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_data_out) {
    // Antenna X, Y, and Z position and reference range to scene center
    auto x_ant_pos_iter = aie::begin(x_ant_pos);
    auto y_ant_pos_iter = aie::begin(y_ant_pos);
    auto z_ant_pos_iter = aie::begin(z_ant_pos);
    auto ref_range_iter = aie::begin(ref_range);

    // Output to backprojection kernel(s)
    auto st_data_out_iter = aie::begin(slowtime_data_out);
    
    *st_data_out_iter++ = *x_ant_pos_iter++;
    *st_data_out_iter++ = *y_ant_pos_iter++;
    *st_data_out_iter++ = *z_ant_pos_iter++;
    *st_data_out_iter++ = *ref_range_iter++;
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

void Backprojection::backprojection_kern(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_data,
                                         input_buffer<TT_DATA, extents<2048>>& __restrict ph_data,
                                         output_async_buffer<TT_DATA, extents<2048>>& __restrict img_out) {
    // Lock img_out buffer
    img_out.acquire();

    // Extract antenna position and range to scene center from slow time data
    auto st_data_vec_iter = aie::begin_vector<ST_ELEMENTS>(slowtime_data);
    auto st_data_vec = *st_data_vec_iter++;

    float x_ant = st_data_vec[0];
    float y_ant = st_data_vec[1];
    float z_ant = st_data_vec[2];
    float r0 = st_data_vec[3];

    // Phase history data in time domain (compressed range lines)
    auto ph_data_iter = aie::begin(ph_data);

    // Image output
    auto img_out_iter = aie::begin_vector<16>(img_out);
    
    // Initialize radar parameters
    init_radar_params(1471301.6, 8192, 9288080400.0);

    // Initialize range grid based on range compressed data
    init_range_grid();
    
    // Initialize target pixel grid
    init_pixel_grid(-5, 30, 32, -5, 5, 32);

    // Initialize pixel segment for this specific kernel
    init_pixel_segment();

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
    alignas(aie::vector_decl_align) float start_x_pxls[16];
    //m_seg.x_pxls = aie::zeros<float,16>();
    //m_seg.y_pxls = aie::zeros<float,16>();
    aie::vector<float,16> x_pxls = aie::zeros<float,16>();
    aie::vector<float,16> y_pxls = aie::zeros<float,16>();
    m_seg.valid_elems = 0;
    if (m_seg.st_x_px_bound) {
        float x_pxl;
        for(int i = 0; i < 16; i++) chess_prepare_for_pipelining {
            x_pxl = m_px_grid.x_st+(m_px_grid.x_res*i);
            if (x_pxl < m_seg.range_en) {
                x_pxls[i] = x_pxl;
                m_seg.valid_elems = i+1;
            }
        }
        for(int i = 0; i < 16; i++) chess_prepare_for_pipelining {
            y_pxls[i] = m_px_grid.y_st+(m_px_grid.y_res*i);
        }
    }
    else if (m_seg.en_x_px_bound || m_seg.all_x_px_bound) {
        // Calculate number of steps for new start position based on x_res (rounds up)
        int steps = (int)((m_seg.range_st-m_px_grid.x_st)/m_px_grid.x_res) + 1;
        float new_x_px_st = m_px_grid.x_st + steps*m_px_grid.x_res;
        float x_pxl;
        for(int i = 0; i < 16; i++) chess_prepare_for_pipelining {
            x_pxl = new_x_px_st+(m_px_grid.x_res*i);
            if (x_pxl < m_seg.range_en) {
                x_pxls[i] = x_pxl;
                m_seg.valid_elems = i+1;
            }
        }
        for(int i = 0; i < 16; i++) chess_prepare_for_pipelining {
            y_pxls[i] = m_px_grid.y_st+(m_px_grid.y_res*i);
        }
    }
    aie::store_v(start_x_pxls, x_pxls);
    
    // Increment x_pxls and y_pxls by this amount every iteration
    float x_px_inc = m_px_grid.x_res*16;
    float y_px_inc = m_px_grid.y_res*16;
    
    if (m_seg.st_x_px_bound || m_seg.en_x_px_bound || m_seg.all_x_px_bound) {
        m_px_grid.display(m_id);
        m_seg.display(m_id);
        m_range_grid.display(m_id);
        printf("%d: x_pxls=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
                x_pxls.get(0),  x_pxls.get(1),  x_pxls.get(2),  x_pxls.get(3), 
                x_pxls.get(4),  x_pxls.get(5),  x_pxls.get(6),  x_pxls.get(7),
                x_pxls.get(8),  x_pxls.get(9),  x_pxls.get(10), x_pxls.get(11), 
                x_pxls.get(12), x_pxls.get(13), x_pxls.get(14), x_pxls.get(15));
        printf("%d: y_pxls=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
                y_pxls.get(0),  y_pxls.get(1),  y_pxls.get(2),  y_pxls.get(3), 
                y_pxls.get(4),  y_pxls.get(5),  y_pxls.get(6),  y_pxls.get(7),
                y_pxls.get(8),  y_pxls.get(9),  y_pxls.get(10), y_pxls.get(11), 
                y_pxls.get(12), y_pxls.get(13), y_pxls.get(14), y_pxls.get(15));
    }
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
    
    aie::vector<TT_DATA,16> low_ph_data;
    aie::vector<TT_DATA,16> high_ph_data;
    if (m_seg.st_x_px_bound || m_seg.en_x_px_bound || m_seg.all_x_px_bound) {
        //printf("ID: %d\n", m_id);
        for(int y_idx=0; y_idx < m_px_grid.y_len/16; y_idx++) chess_prepare_for_pipelining {
            for(int x_idx=0; x_idx < m_px_grid.x_len/16; x_idx++) chess_prepare_for_pipelining {

                // Calculate differential range for pixel segments
                auto x = aie::sub(x_ant, x_pxls);
                auto x_diff_sq = aie::mul_square(x).to_vector<float>(0);

                auto y = aie::sub(y_ant, y_pxls);
                auto y_diff_sq = aie::mul_square(y).to_vector<float>(0);

                auto xy_add = aie::add(x_diff_sq, y_diff_sq);
                auto xyz_add = aie::add(xy_add, z_diff_sq);

                auto xyz_sqrt = aie::sqrt(xyz_add);

                auto differ_range = aie::sub(xyz_sqrt, r0);

                // Calculate phase correction for image
                auto ph_corr_angle = aie::mul(m_radar_params.ph_corr_coef, differ_range).to_vector<float>(0);
                auto ph_corr = aie::sincos_complex(ph_corr_angle);
                
                // Determine which pixels fall within the swath range
                //auto right_mask = aie::ge(differ_range, min_range_bins);
                //auto left_mask = aie::le(differ_range, max_range_bins);
                //auto differ_range_mask = right_mask & left_mask;
                
                //TODO: WRONG
                auto px_diff = aie::sub(m_range_grid.range_width, x_pxls);
                printf("%d: px_diff=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
                        px_diff.get(0),  px_diff.get(1),  px_diff.get(2),  px_diff.get(3), 
                        px_diff.get(4),  px_diff.get(5),  px_diff.get(6),  px_diff.get(7),
                        px_diff.get(8),  px_diff.get(9),  px_diff.get(10), px_diff.get(11), 
                        px_diff.get(12), px_diff.get(13), px_diff.get(14), px_diff.get(15));

                //aie::saturation_mode current_sat=aie::get_saturation();
                //printf("CURRENT_SATURATION = %d\n", current_sat);

                auto low_range_bins_float = aie::mul(px_diff, m_range_grid.inv_range_res).to_vector<float>(0);
                printf("%d: low_range_bins_float=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
                        low_range_bins_float.get(0),  low_range_bins_float.get(1),  low_range_bins_float.get(2),  low_range_bins_float.get(3), 
                        low_range_bins_float.get(4),  low_range_bins_float.get(5),  low_range_bins_float.get(6),  low_range_bins_float.get(7),
                        low_range_bins_float.get(8),  low_range_bins_float.get(9),  low_range_bins_float.get(10), low_range_bins_float.get(11), 
                        low_range_bins_float.get(12), low_range_bins_float.get(13), low_range_bins_float.get(14), low_range_bins_float.get(15));

                // TODO: This is actually rounding correctly...but I need it to floor round...
                aie::vector<int32,16> low_range_bins = aie::to_fixed<int32>(low_range_bins_float);
                printf("%d: low_range_bins=%f %f %f %f %f %f %f %f %f %f %f %f %f %f %f %f\n", m_id, 
                        low_range_bins.get(0),  low_range_bins.get(1),  low_range_bins.get(2),  low_range_bins.get(3), 
                        low_range_bins.get(4),  low_range_bins.get(5),  low_range_bins.get(6),  low_range_bins.get(7),
                        low_range_bins.get(8),  low_range_bins.get(9),  low_range_bins.get(10), low_range_bins.get(11), 
                        low_range_bins.get(12), low_range_bins.get(13), low_range_bins.get(14), low_range_bins.get(15));

                //low_range_bins = aie::sub(low_range_bins, fft_seg);
                //auto high_range_bins = aie::add(low_range_bins, 1);

                
                //for (int i=0; i<16; i++) {
                //    low_ph_data[i] = ph_data_iter[low_range_bins[i]];
                //    high_ph_data[i] = ph_data_iter[low_range_bins[i]+1];
                //    printf("ph_corr[%d] = {%f, %f}\n", i, ((cfloat)ph_corr[i]).real, ((cfloat)ph_corr[i]).imag);
                //    printf("low_ph_data[%d] = {%f, %f}\n", i, ((cfloat)low_ph_data[i]).real, ((cfloat)low_ph_data[i]).imag);
                //    printf("high_ph_data[%d] = {%f, %f}\n", i, ((cfloat)high_ph_data[i]).real, ((cfloat)high_ph_data[i]).imag);
                //}

                //y0+(x-x0)*((y1 - y0)/dr);
                // Interpolation...HOW TO DO THIS
                //y0+(x-x0)*((y1 - y0)/(x1-x0));
                
                //if (differ_range_mask.full()) {
                //    auto img_vec = aie::mul(differ_range, ph_corr).to_vector<cfloat>(0);
                //    aie::store_v(img+(y_idx*x_px_grid_len)+(x_idx*16), img_vec);
                //}

                //// Increment range bins to next segment
                //range_bins = aie::add(dr, range_bins);
                
                // Increment X pixel segment
                x_pxls = aie::add(x_pxls, x_px_inc);
            }

            // Reset X pixel segment
            x_pxls = aie::load_v<16>(start_x_pxls);

            // Increment Y pixel segment
            y_pxls = aie::add(y_pxls, y_px_inc);
        }
    }

    //img_out.release();
}

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
