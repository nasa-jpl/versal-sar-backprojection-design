// By: Austin Owens
// Date: 6/6/2024
// Desc: Header for synchronizing variables across projects
#ifndef COMMON_H
#define COMMON_H

// Macro to convert a token to a string
#define STRINGIFY(x) #x
#define TO_STRING(x) STRINGIFY(x)

// Number of pulses to process. This will also be the number of rows of the
// output image. Because the work to be processed is divided equally across all
// AI cores, this number must be selected carefully so there is no leftover
// work to be processed. As a result, this equation (i.e. pixels to be
// processed per img solver) must be equal to a whole number:
// (RC_SAMPLES * PULSES) / (IMG_SOLVERS_PER_SWITCH * AIE_SWITCHES)
#define PULSES 602

// Number of range compression samples. This will also be the number of columns
// of the output image. Test data only generated for below available options
// AVAILABLE OPTIONS: 512, 256, 128, 64 
#define RC_SAMPLES 512

// Number of azimuth samples. This will also be the number of rows of the
// output image. TODO: It would be nice to add AZ_SAMPLES so that PULSES and
// AZ_SAMPLES can be differentiated when the AIE engine is performing the
// cumulative summation over the AZ_SAMPLES. 
//#define AZ_SAMPLES PULSES

// Number of AIE switches to use
#define AIE_SWITCHES 7

// Number of image reconstruction solvers per switch. It must be a power of 2
// and the max number of image solvers on a switch is 32
#define IMG_SOLVERS_PER_SWITCH 32

// Number of image reconstruction solvers (must be pwr of 2)
#define IMG_SOLVERS (AIE_SWITCHES*IMG_SOLVERS_PER_SWITCH)

// Number of broadcasted elements to other AI kernels:
// X position of antenna 
// Y position of antenna 
// Z position of antenna 
// Range to center of scene (ref_range)
#define BC_ELEMENTS 4

// Constants
static constexpr float PI = 3.1415926535898;
static constexpr float TWO_PI = 6.2831853071796;
static constexpr float INV_TWO_PI = 0.1591549430919; // Used in AIE code only

// Radar parameters
static constexpr float C = 299792458.0;
static constexpr float MIN_FREQ = 9288080400.0; // Used in AIE code only
static constexpr float RANGE_FREQ_STEP = 1471301.6;
static constexpr float RANGE_WIDTH = C/(2.0*RANGE_FREQ_STEP);
static constexpr float RANGE_RES = RANGE_WIDTH/RC_SAMPLES;
static constexpr float INV_RANGE_RES = 1.0/RANGE_RES; // Used in AIE code only
static constexpr int HALF_RANGE_SAMPLES = RC_SAMPLES/2;

#endif // COMMON_H
