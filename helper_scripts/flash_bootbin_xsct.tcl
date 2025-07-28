# By: Austin Owens
# Date: 9/3/2024
# Desc: Helper tcl script to flash BOOT.BIN into versal over JTAG

set design_build_path [file dirname [file dirname [file normalize [info script]]]]

# Connect
connect

# Select a target
targets -set -nocase -filter {name =~ "Versal*"}

# System Reset
rst

# Program the file
device program $design_build_path/build/hw/package/BOOT.BIN
