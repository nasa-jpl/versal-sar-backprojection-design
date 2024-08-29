// By: Austin Owens
// Date: 6/6/2024
// Desc: Header for synchronizing variables across projects
#ifndef COMMON_H
#define COMMON_H

// Macro to convert a token to a string
#define STRINGIFY(x) #x
#define TO_STRING(x) STRINGIFY(x)

// FFT VARIABLES

// Influences TT_DATA and TT_TWIDDLE. It is the Data type of input into AIE. 
// cint16 (0) or cfloat (1)
#define DATA_TYPE 1

// Is an unsigned integer which describes the number of samples in the transform.
// This must be 2^N where N is an integer in the range 4 to 16 inclusive.
// When TP_DYN_PT_SIZE is set, TP_POINT_SIZE describes the maximum point size possible.
#define TP_POINT_SIZE 1024
#define MAT_ROWS TP_POINT_SIZE
#define MAT_COLS TP_POINT_SIZE
#define TMPL_MAT_ROWS TP_POINT_SIZE/8
#define TMPL_MAT_COLS TP_POINT_SIZE/8

// Selects whether the transform to perform is an FFT (1) or IFFT (0).
#define TP_FFT 1
#define TP_IFFT 0

// Selects the power of 2 to scale the result by prior to output.
#define TP_FFT_SHIFT 0

// Selects whether (1) or not (0) to use run-time point size determination.
// When set, each window of data must be preceeded, in the window, by a 256 bit header.
// This header is 8 samples when TT_DATA is cint16 and 4 samples otherwise.
// The real part of the first sample indicates the forward (1) or inverse (0) transform.
// The real part of the second sample indicates the Radix2 power of the point size.
//
// e.g. for a 512 point size, this field would hold 9, as 2^9 = 512. The second least 
// significant byte 8 bits of this field describe the Radix 2 power of the following
// frame. e.g. for a 512 point size, this field would hold 9, as 2^9 = 512.
//
// Any value below 4 or greater than log2(TP_POINT_SIZE) is considered illegal.
// The output window will also be preceeded by a 256 bit vector which is a copy of the 
// input vector, but for the real part of the top sample, which is 0 to indicate a legal 
// frame or 1 to indicate an illegal frame. When TP_PARALLEL_POWER is greater than 0, the 
// header must be applied before each window of data for every port of the design and will 
// appears before each window of data on the output ports. 
//
// Note that the minimum point size of 16 applies to each lane when in parallel mode, 
// so a configuration of point size 256 with TP_PARALLEL_POWER = 2 will have 4 lanes 
// each with a minimum of 16 so the minimum legal point size here is 64.
#define TP_DYN_PT_SIZE 0

// Is an unsigned integer which describes the number of samples to be processed in 
// each call to the function. When TP_DYN_PT_SIZE is set to 1 the actual window size
// will be larger than TP_FFT_WINDOW_VSIZE because the header is not included in 
// TP_FFT_WINDOW_VSIZE. By default, TP_FFT_WINDOW_VSIZE is set to match TP_POINT_SIZE.
// TP_FFT_WINDOW_VSIZE may be set to be an integer multiple of the TP_POINT_SIZE, in 
// which case multiple FFT iterations will be performed on a given input window, 
// resulting in multiple iterations of output samples, reducing the numer of times the 
// kernel needs to be triggered to process a given number of input data samples. As a 
// result, the overheads inferred during kernel triggering are reduced and overall 
// performance is increased.
#define TP_FFT_WINDOW_VSIZE TP_POINT_SIZE

#if DATA_TYPE == 0
    // Used for the aiesim to validate data at various stages of the pipeline
    #define FFT_SAMPLE_DATA 1
    #define HP_SAMPLE_DATA 2

    // Describes the type of individual data samples input to and output from the
    // transform function. This is a typename and must be one of the following:
    // int16, cint16, int32, cint32, float, cfloat.
    #define TT_DATA cint16

    // Describes the type of twiddle factors of the transform. It must be one of the
    // following: cint16, cint32, cfloat and must also satisfy the following rules:
    // - 32 bit types are only supported when TT_DATA is also a 32 bit type
    // - TT_TWIDDLE must be an integer type if TT_DATA is an integer type
    // - TT_TWIDDLE must be cfloat type if TT_DATA is a float type
    #define TT_TWIDDLE cint16

    // Kernel I/O window buff size in bytes
    #define FFT_WINDOW_BUFF_SIZE (TP_FFT_WINDOW_VSIZE * 4)

    // Selects the number of kernels the FFT will be divided over in series to improve 
    // throughput
    #define TP_FFT_CASC_LEN 1
   
#elif DATA_TYPE == 1
    // Used for the aiesim to validate data at various stages of the pipeline
    #define FFT_SAMPLE_DATA 1.5
    #define HP_SAMPLE_DATA 4.5

    // Describes the type of individual data samples input to and output from the
    // transform function. This is a typename and must be one of the following:
    // int16, cint16, int32, cint32, float, cfloat.
    #define TT_DATA cfloat

    // Describes the type of twiddle factors of the transform. It must be one of the
    // following: cint16, cint32, cfloat and must also satisfy the following rules:
    // - 32 bit types are only supported when TT_DATA is also a 32 bit type
    // - TT_TWIDDLE must be an integer type if TT_DATA is an integer type
    // - TT_TWIDDLE must be cfloat type if TT_DATA is a float type
    #define TT_TWIDDLE cfloat
    
    // Kernel I/O window buff size in bytes
    #define FFT_WINDOW_BUFF_SIZE (TP_FFT_WINDOW_VSIZE * 8)

    // Selects the number of kernels the FFT will be divided over in series to improve 
    // throughput
    #define TP_FFT_CASC_LEN 4
   
#endif

// HADAMARD PRODUCT VARIABLES

// Describes the type of individual data samples input to the function. This is a 
// typename and must be one of the following:
// int16, int32, cint16, cint32, float, cfloat.
#define TT_DATA_A TT_DATA

// Describes the type of individual data samples input to the function. This is a 
// typename and must be one of the following:
// int16, int32, cint16, cint32, float, cfloat.
#define TT_DATA_B TT_DATA

// Describes the number of samples in the vectors A and B.
#define TP_DIM TP_POINT_SIZE

// Describes the number of vectors to be processed in each call to this function.
#define TP_NUM_FRAMES 4

// Describes power of 2 shift down applied to the accumulation of product terms
// before each output. TP_SHIFT must be in the range 0 to 61.
#define TP_HP_SHIFT 0

// Described whether to use streams (1) or windows (0).
#define TP_API 0

// Describes the number of kernels to use in parallel.
#define TP_SSR 2

// Describes the selection of rounding to be applied during the shift down stage of
// processing. Although, TP_RND accepts unsigned integer values descriptive macros are
// recommended where:
// - rnd_floor = Truncate LSB, always round down (towards negative infinity).
// - rnd_ceil = Always round up (towards positive infinity).
// - rnd_sym_floor = Truncate LSB, always round towards 0.
// - rnd_sym_ceil = Always round up towards infinity.
// - rnd_pos_inf = Round halfway towards positive infinity.
// - rnd_neg_inf = Round halfway towards negative infinity.
// - rnd_sym_inf = Round halfway towards infinity (away from zero).
// - rnd_sym_zero = Round halfway towards zero (away from infinity).
// - rnd_conv_even = Round halfway towards nearest even number.
// - rnd_conv_odd = Round halfway towards nearest odd number.
// No rounding is performed on ceil or floor mode variants. Other modes round to the
// nearest integer. They differ only in how they round for values of 0.5.
#define TP_RND 6

// Describes the selection of saturation to be applied during the shift down stage 
// of processing. TP_SAT accepts unsigned integer values, where:
//
// 0: none = No saturation is performed and the value is truncated on the MSB side.
//
// 1: saturate = Default. Saturation rounds an n-bit signed value in the range 
// [- ( 2^(n-1) ) : +2^(n-1) - 1 ].
//
// 3: symmetric = Controls symmetric saturation. Symmetric saturation rounds an n-bit 
// signed value in the range [- ( 2^(n-1) -1 ) : +2^(n-1) - 1 ].
#define TP_SAT 1

#endif // COMMON_H
