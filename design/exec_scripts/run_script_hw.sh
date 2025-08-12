#!/bin/bash
export XILINX_XRT=/usr

# Check if exactly one argument is provided
if [ "$#" -ne 4 ]; then
	echo "Usage: $0 <slowtime dataset csv file> <range compressed dataset csv file> <img out csv file> <iteration>"
	exit 1
fi

# Get the data type and matrix dimension from the arguments
st_csv_file=$1
rc_csv_file=$2
img_out_csv_file=$3
iter=$4
./sar_backproject.elf a.xclbin $st_csv_file $rc_csv_file $img_out_csv_file $iter
