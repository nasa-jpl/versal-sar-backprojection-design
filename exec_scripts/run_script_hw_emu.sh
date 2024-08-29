#!/bin/bash
export XILINX_XRT=/usr
export XCL_EMULATION_MODE=hw_emu

# Check if exactly two arguments are provided
if [ "$#" -ne 2 ]; then
	echo "Usage: $0 <data_type> <mat_dim>"
	exit 1
fi

# Get the data type and matrix dimension from the arguments
data_type=$1
mat_dim=$2

# Validate the data type
if [[ "$data_type" != "cint16" && "$data_type" != "cfloat" ]]; then
	echo "Error: data_type must be either 'cint16' or 'cfloat'"
	exit 1
fi

# Construct the file path with the provided data type and matrix dimension
file_path="test_data/${data_type}_simple_mat/mat_${mat_dim}x${mat_dim}.txt"

# Run the command with the constructed file path
./img_fft_cross_corr.elf a.xclbin "$file_path" "$file_path"
