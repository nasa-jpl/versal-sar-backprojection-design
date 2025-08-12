# Run this with the following commnad: 'vitis-run --tcl run_dma_pkt_router_tb.tcl'

# Create a project
open_project -reset dma_pkt_router_testbench

# Change working directory to where this tcl script is located
set script_dir [file dirname [info script]]

# Add design files
add_files $script_dir/../dma_pkt_router.cpp

# Add test bench & files
add_files -tb $script_dir/dma_pkt_router_tb.cpp

# Set the top-level function
set_top dma_pkt_router

# ########################################################

# Create a solution
open_solution -reset solution1 -flow_target vitis

# Define technology and clock rate
set_part  {xcvc1902-vsva2197-2MP-e-S}
create_clock -period "312.5MHz"

# Set variable to select which steps to execute
set hls_exec 2

csim_design

# Set any optimization directives
# End of directives

if {$hls_exec == 1} {
	# Run Synthesis and Exit
	csynth_design
	
} elseif {$hls_exec == 2} {
	# Run Synthesis, RTL Simulation and Exit
	csynth_design
	
	cosim_design
} elseif {$hls_exec == 3} { 
	# Run Synthesis, RTL Simulation, RTL implementation and Exit
	csynth_design
	
	cosim_design
	export_design -format ip_catalog
} else {
	# Default is to exit after setup
	csynth_design
}

exit


