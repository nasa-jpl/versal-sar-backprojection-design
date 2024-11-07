// By: Austin Owens
// Date: 6/25/2024
// Desc: Kernel header

#pragma once

#include <adf.h>
#include "aie_api/aie.hpp"
#include "../common.h"

using namespace adf;

void cplx_conj_kern(input_buffer<TT_DATA, extents<2048>>& __restrict in, 
                    output_buffer<TT_DATA, extents<2048>>& __restrict out);

void slowtime_splicer_kern(input_buffer<float, extents<1>>& __restrict x_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict y_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict z_ant_pos_in,
                           input_buffer<float, extents<1>>& __restrict ref_range_in,
                           output_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_out);

class Backprojection
{
    public:
        static constexpr float PI = 3.1415926535898;
        static constexpr float TWO_PI = 6.2831853071796;
        static constexpr float INV_TWO_PI = 0.15915494309189;
        static constexpr float C = 299792458.0;
        static constexpr int ACCUM_PULSES = 2;

        Backprojection(int id);
        void backprojection_kern(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_in,
                                 input_circular_buffer<TT_DATA, extents<2048>>& __restrict rc_in,
                                 input_buffer<TT_DATA, extents<2048>>& __restrict xy_px_in,
                                 output_async_buffer<TT_DATA, extents<2048>>& __restrict img_out);

        void generate_pixel_grid(float x_st, float x_en, float y_st, float y_en, float x_grid_len, float y_grid_len);
        static void registerKernelClass()
        { 
            REGISTER_FUNCTION(Backprojection::backprojection_kern);
        }

    private:
        uint32 m_id;
        uint32 m_iter;
        alignas(aie::vector_decl_align) TT_DATA m_img[2048];

        void init_radar_params(float range_freq_step, int range_samples, float min_freq);
        void init_range_grid();
        void init_pixel_grid(float x_st, float x_en, int x_len, float y_st, float y_en, int y_len);
        void init_pixel_segment();
        void print_float(const char *str, aie::vector<float,16> data);
        void print_int32(const char *str, aie::vector<int32,16> data);
        
        // Radar parameters
        struct RadarParams {
            float range_freq_step; // Range frequency step size
            int range_samples;     // Number of samples per range line
            float min_freq;        // Minimum frequency in range line
            float ph_corr_coef;    // Phase correction coefficient
            void display(int id) const {
                printf("\nRADAR PARAMS (ID=%d)\nrange_freq_step: %f\nrange_samples: %d\n" \
                        "min_freq: %f\nph_corr_coef: %f\n", id, range_freq_step, range_samples, 
                                                            min_freq, ph_corr_coef);
            }
        } m_radar_params;

        // Range grid based on signal time delay
        struct RangeGrid {
            float range_width;     // Max width scene size of image (in range direction)
            float range_res;       // Range step size
            float inv_range_res;   // Inverse range step size
            float seg_offset;      // Segment offset of samples
            void display(int id) const {
                printf("\nRANGE GRID DATA (ID=%d)\nrange_width: %f\nrange_res: %f\n", id, range_width, range_res);
            }
        } m_range_grid;
        
        //// Pixel grid representing the desired target
        //struct PixelGrid {
        //    float x_st;  // Desired pixel start target in x direction
        //    float x_en;  // Desired pixel end target in x direction
        //    int x_len;   // Desired pixel length in x direction
        //    float x_res; // Pixel x resolution
        //    float y_st;  // Desired pixel start target in y direction
        //    float y_en;  // Desired pixel end target in y direction
        //    int y_len;   // Desired pixel length in y direction
        //    float y_res; // Pixel y resolution
        //    void display(int id) const {
        //        printf("\nPIXEL GRID DATA (ID=%d):\n" \
        //               "x_st: %f | x_en: %f | x_len: %d | x_res: %f\n" \
        //               "y_st: %f | y_en: %f | y_len: %d | y_res: %f\n", id, x_st, x_en, x_len, x_res, 
        //                                                                y_st, y_en, y_len, y_res);
        //    }
        //} m_px_grid;

        //// Range/pixel segment values based on instantiated backprojection kernel
        //struct RangePixelSegment {
        //    float range_st;        // The starting range based on current segment
        //    float range_en;        // The ending range based on current segment
        //    bool st_x_px_bound;    // Bool if starting X pixel is bounded in range_st and range_en
        //    bool en_x_px_bound;    // Bool if ending X pixel is bounded in range_st and range_en
        //    bool all_x_px_bound;   // Bool if all X pixels are bounded in range_st and range_en
        //    int valid_elems;       // Valid elemnts in x_pxls vector
        //    void display(int id) const {
        //        printf("\nRANGE PIXEL SEGMENT DATA (ID=%d)\nrange_st: %f | range_en: %f\n" \
        //                "st_x_px_bound: %d | en_x_px_bound: %d | all_x_px_bound: %d\nvalid_elems: %d\n", 
        //                id, range_st, range_en, st_x_px_bound, en_x_px_bound, all_x_px_bound, valid_elems);
        //    }
        //} m_seg;

};

//void backprojection_kern(input_buffer<float, extents<ST_ELEMENTS>>& __restrict slowtime_data,
//                         input_buffer<TT_DATA, extents<2048>>& __restrict ph_data,
//                         output_buffer<TT_DATA, extents<2048>>& __restrict img_out,
//                         const int instance);

