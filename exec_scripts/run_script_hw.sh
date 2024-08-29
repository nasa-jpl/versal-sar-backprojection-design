#!/bin/bash
export XILINX_XRT=/usr

## Check if exactly two arguments are provided
#if [ "$#" -ne 2 ]; then
#	echo "Usage: $0 <data_type> <mat_dim>"
#	exit 1
#fi
#
## Get the data type and matrix dimension from the arguments
#data_type=$1
#mat_dim=$2
#
## Validate the data type
#if [[ "$data_type" != "cint16" && "$data_type" != "cfloat" ]]; then
#	echo "Error: data_type must be either 'cint16' or 'cfloat'"
#	exit 1
#fi
#
## Construct the file path with the provided data type and matrix dimension
#file_path="test_data/${data_type}_simple_mat/mat_${mat_dim}x${mat_dim}.txt"
#
## Run the command with the constructed file path
#./img_fft_cross_corr_xrt.elf a.xclbin "$file_path" "$file_path"

# Check if exactly one argument is provided
if [ "$#" -ne 1 ]; then
	echo "Usage: $0 <iterations>"
	exit 1
fi

# Get the data type and matrix dimension from the arguments
iter=$1
./img_fft_cross_corr.elf a.xclbin ./test_data/test01_map.hex "$iter"
