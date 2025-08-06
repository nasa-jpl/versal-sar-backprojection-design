# By: Austin Owens
# Date: 5/21/2024
# Desc: Makefile for building various elements of the platform for image FFT cross-correlation
#
# Targets: package, run, pl, plsim_router, aie, host, aiesim, aiesim_profile, aiesim_xpe, metrics, clean
# 
# Examples:
# make
# make aie
# make aiesim
# make host
# make package
# make TARGET=sw_emu run
# make TARGET=hw_emu run


##### VARIABLES #####

# User Defined Variables
TARGET ?= hw

# Global Vars
XSA = sar_backproject_${TARGET}.xsa
HOST_EXE = sar_backproject.elf
PROJECT_DIR = $(patsubst %/,%, $(dir $(realpath $(firstword $(MAKEFILE_LIST)))))
EMU_LAUNCH_FILE = launch_${TARGET}.sh
ifeq ($(TARGET), sw_emu)
    AIE_TARGET = x86sim
    PL_TARGET = x86
else
    AIE_TARGET = hw
    PL_TARGET = hw
endif

# Build Directories
BUILD_DIR = build/${TARGET}
AIE_BUILD_DIR = ${BUILD_DIR}/aie/${AIE_TARGET}
PL_BUILD_DIR = ${BUILD_DIR}/pl/${TARGET}
AIESIM_BUILD_DIR = ${BUILD_DIR}/aiesim
PLSIM_BUILD_DIR = ${BUILD_DIR}/plsim
HOST_BUILD_DIR = ${BUILD_DIR}/host
XSA_BUILD_DIR = ${BUILD_DIR}/xsa
PACKAGE_BUILD_DIR = ${BUILD_DIR}/package
METRICS_BUILD_DIR = build/hw/metrics

##### CHECKS #####

# Check valid TARGET value
VALID_TARGETS := hw hw_emu sw_emu
ifeq ($(filter $(TARGET),$(VALID_TARGETS)),)
    $(error Invalid TARGET specified: $(TARGET). Valid options are 'hw', 'hw_emu', or 'sw_emu')
endif

# Check if PLATFORM is set from the environment
ifndef PLATFORM
    $(error PLATFORM is not set. Please source the env_setup.sh script.)
endif

# Check for incompatible TARGET=sw_emu and 'aiesim' combination
ifneq (,$(filter aiesim,$(MAKECMDGOALS)))
ifeq ($(TARGET), sw_emu)
    $(error The 'aiesim' target cannot be used with 'sw_emu'. \
The AI Engine simulator (aiesimulator) models the timing and resources of \
the AI Engine array and requires the application to be compiled with \
hardware (hw) or hardware emulation (hw_emu) targets to provide \
cycle-accurate performance analysis. The 'sw_emu' target uses the x86 \
simulator, which provides a fast functional simulation without timing \
accuracy. Please use 'hw' or 'hw_emu' for accurate AI Engine simulation)
endif
endif

# Check for incompatible TARGET=hw and 'run' combination
ifneq (,$(filter run,$(MAKECMDGOALS)))
ifeq ($(TARGET), hw)
    $(error The 'run' target can only be used with a TARGET=hw_emu or TARGET=sw_emu)
endif
endif

.PHONY: package run pl plsim_router aie aiesim aiesim_profile host aiesim_xpe metrics clean


##### DIRECT TARGETS #####

# Packaging target for linking everything together (only package the DTB with hw since
# custom DTB's seems to cause the kernel to panic and crash when QEMU is used)
ifeq ($(TARGET),hw)
    DTB_OPTION := --package.dtb ${DTB}
else
    DTB_OPTION :=
endif
package: ${PL_BUILD_DIR}/dma_pkt_router.xo ${AIE_BUILD_DIR}/libadf.a ${HOST_BUILD_DIR}/${HOST_EXE} ${XSA_BUILD_DIR}/${XSA} 
	@RC_SAMPLES=$$(grep '^#define RC_SAMPLES' ${PROJECT_DIR}/design/common.h | awk '{print $$3}'); \
	mkdir -p ${PACKAGE_BUILD_DIR}; \
	cd ${PACKAGE_BUILD_DIR}; \
	v++ -p -t ${TARGET} -f ${PLATFORM} \
		$(DTB_OPTION) \
		--package.bl31_elf ${BL31_ELF} \
		--package.uboot ${UBOOT} \
		--package.kernel_image ${IMAGE} \
		--package.rootfs ${ROOTFS} \
		--package.boot_mode=sd \
		--package.image_format=ext4 \
		--package.defer_aie_run \
		--package.sd_file ${PROJECT_DIR}/design/exec_scripts/run_script_${TARGET}.sh \
		--package.sd_file ${PROJECT_DIR}/design/profiling_cfgs/xrt.ini \
		--package.sd_file ${PROJECT_DIR}/${HOST_BUILD_DIR}/${HOST_EXE} \
		--package.sd_file ${PROJECT_DIR}/design/test_data/gotcha_slowtime_pass1_360deg_HH.csv \
		--package.sd_file ${PROJECT_DIR}/design/test_data/gotcha_$${RC_SAMPLES}-out-of-424-rc-samples_pass1_360deg_HH.csv \
		${PROJECT_DIR}/${XSA_BUILD_DIR}/${XSA} ${PROJECT_DIR}/${AIE_BUILD_DIR}/libadf.a
	@echo ""
	@echo "Packaging Complete..."
	@echo "####################################"
	@echo ""

# Run target for running emulation. Influenced by TARGET being sw_emu or hw_emu.
run: package
	${PACKAGE_BUILD_DIR}/${EMU_LAUNCH_FILE}
	@echo ""
	@echo "Emulation Finished..."
	@echo "####################################"
	@echo ""

# PL target for building all HLS Programmable Logic (PL) kernels
pl: ${PL_BUILD_DIR}/dma_pkt_router.xo

# PL sim target for building all artifacts needed for running the stride
# controller PL sim.
#plsim_stride:
#	mkdir -p ${PLSIM_BUILD_DIR}/plsimulator_output; \
#	cd ${PLSIM_BUILD_DIR}; \
#	vitis-run --tcl ${PROJECT_DIR}/design/pl/tb/run_dma_stride_controller_tb.tcl | tee plsim_stride.log
#	@echo ""
#	@echo "DMA Stride Controller PL Simulation, Complete..."
#	@echo "####################################"
#	@echo ""

# PL sim target for building all artifacts needed for running the packet router
# PL sim. If the CSV file doesn't exist, run aiesim to generate it.
plsim_router:
	if [ ! -f ${AIESIM_BUILD_DIR}/aiesimulator_output/aie_to_plio_switch_0_0.csv ]; then \
	    $(MAKE) aiesim; \
	fi
	mkdir -p ${PLSIM_BUILD_DIR}/plsimulator_output; \
	cd ${PLSIM_BUILD_DIR}; \
	vitis-run --tcl ${PROJECT_DIR}/design/pl/tb/run_dma_pkt_router_tb.tcl | tee plsim_router.log
	@echo ""
	@echo "DMA Packet Router PL Simulation, Complete..."
	@echo "####################################"
	@echo ""


# AIE target for building all AIE related apps
aie: ${AIE_BUILD_DIR}/libadf.a

# Host target for building all host related apps
host: ${HOST_BUILD_DIR}/${HOST_EXE}

# AIE simulation target
aiesim: ${AIE_BUILD_DIR}/libadf.a
	mkdir -p ${AIESIM_BUILD_DIR}; \
	cd ${AIESIM_BUILD_DIR}; \
	aiesimulator --pkg-dir=${PROJECT_DIR}/${BUILD_DIR}/Work \
		--input-dir ${PROJECT_DIR}/${PLSIM_BUILD_DIR}/plsimulator_output 2>&1 | tee aiesim.log
	@echo ""
	@echo "AIE Simulation, Without Profiling, Complete..."
	@echo "####################################"
	@echo ""

# AIE simulation target that enables profiling. This allows prints to appear in
# console in AIE kernel code and also generates profile information.
aiesim_profile: ${AIE_BUILD_DIR}/libadf.a
	mkdir -p ${AIESIM_BUILD_DIR}; \
	cd ${AIESIM_BUILD_DIR}; \
	aiesimulator --profile --pkg-dir=${PROJECT_DIR}/${BUILD_DIR}/Work \
		--dump-vcd aie \
		--input-dir ${PROJECT_DIR}/${PLSIM_BUILD_DIR}/plsimulator_output 2>&1 | tee aiesim.log
	@echo ""
	@echo "AIE Simulation, With Profiling, Complete..."
	@echo "####################################"
	@echo ""

# Generates power metrics for AIE using the aie.vcd generated by the aiesimulator. The xpe
# file can be used in addition to the xpe file generated from the 'metrics' target to 
# generate a more comprehensive power estimation. You can use the PDM GUI to import these 
# files to look at these various power-related metrics.
aiesim_xpe: ${AIE_BUILD_DIR}/libadf.a ${AIESIM_BUILD_DIR}/aie.vcd
	mkdir -p ${METRICS_BUILD_DIR}; \
	cd ${METRICS_BUILD_DIR}; \
	vcdanalyze --vcd ${PROJECT_DIR}/${AIESIM_BUILD_DIR}/aie.vcd \
		--pkg-dir=${PROJECT_DIR}/${BUILD_DIR}/Work \
		--xpe
	@echo ""
	@echo "AIE XPE Generation Complete..."
	@echo "####################################"
	@echo ""

# Generate utilization and power metrics from entire project. TARGET must be HW.
metrics: ${XSA_BUILD_DIR}/${XSA}
	mkdir -p $(METRICS_BUILD_DIR); \
	cd $(METRICS_BUILD_DIR); \
	vivado -mode batch -source ${PROJECT_DIR}/design/vivado_metrics_scripts/report_metrics.tcl \
		${PROJECT_DIR}/$(XSA_BUILD_DIR)/_x/link/vivado/vpl/prj/prj.xpr
	@echo ""
	@echo "Vivado Utilization/Power Report Generation Complete..."
	@echo "####################################"
	@echo ""

##### INDIRECT TARGETS #####

# Building DMA Stride Controller PL kernel
#${PL_BUILD_DIR}/dma_stride_controller.xo: design/pl/dma_stride_controller.cpp design/pl/dma_stride_controller.h design/pl/stride_controller_config.cfg ${PROJECT_DIR}/design/common.h
#	mkdir -p ${PL_BUILD_DIR}; \
#	cd ${PL_BUILD_DIR}; \
#	v++ -c --mode hls --platform=${PLATFORM} -t ${PL_TARGET} \
#		--work_dir=${PROJECT_DIR}/${BUILD_DIR}/Work \
#		--config ${PROJECT_DIR}/design/pl/stride_controller_config.cfg; \
#	mv ${PROJECT_DIR}/${BUILD_DIR}/Work/dma_stride_controller.xo ./
#	@echo ""
#	@echo "DMA Stride Controller PL Kernel Built..."
#	@echo "####################################"
#	@echo ""

# Building DMA Packet Router PL kernel
${PL_BUILD_DIR}/dma_pkt_router.xo: design/pl/dma_pkt_router.cpp design/pl/dma_pkt_router.h design/pl/pkt_router_config.cfg ${PROJECT_DIR}/design/common.h
	@AIE_SWITCHES=$$(grep '^#define AIE_SWITCHES' ${PROJECT_DIR}/design/common.h | awk '{print $$3}'); \
	echo "# THIS FILE IS AUTO-GENERATED FROM THE MAKEFILE" > ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
	echo -e "\n[clock]" >> ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
	echo -e "default_freqhz=312500000\n" >> ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
	echo "[connectivity]" >> ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
	echo -e "\n### DMA PACKET ROUTER CONTROLLER ###" >> ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
    routers=$$(printf "dma_pkt_router_%s," $$(seq 0 $$((AIE_SWITCHES - 1)))); \
	routers=$${routers%,}; \
	echo "nk=dma_pkt_router:$$AIE_SWITCHES:$$routers" >> ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
	echo -e "\n# Connect AIE graph's plio_pkt_rtr_out_0_# to PL kernel's pl_stream_in" >> ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
	for i in $$(seq 0 $$((AIE_SWITCHES - 1))); do \
		echo "stream_connect=ai_engine_0.plio_pkt_rtr_out_0_$$i:dma_pkt_router_$$i.pl_stream_in" >> ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
	done; \
	echo -e "\n# System port connection linking dma_pkt_router_0 instance to mem resource" >> ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
	for i in $$(seq 0 $$((AIE_SWITCHES - 1))); do \
		echo "sp=dma_pkt_router_$$i.ddr_mem:DDR" >> ${PROJECT_DIR}/design/system_cfgs/system.cfg; \
	done; \
	mkdir -p ${PL_BUILD_DIR}; \
	cd ${PL_BUILD_DIR}; \
	v++ -c --mode hls --platform=${PLATFORM} -t ${PL_TARGET} \
		--work_dir=${PROJECT_DIR}/${BUILD_DIR}/Work \
		--config ${PROJECT_DIR}/design/pl/pkt_router_config.cfg; \
	mv ${PROJECT_DIR}/${BUILD_DIR}/Work/dma_pkt_router.xo ./
	@echo ""
	@echo "DMA Packet Router PL Kernel Built..."
	@echo "####################################"
	@echo ""

# Building the ADF library
${AIE_BUILD_DIR}/libadf.a: design/aie/* ${PROJECT_DIR}/design/common.h design/aie/aiecompiler.cfg
	mkdir -p ${AIE_BUILD_DIR}; \
	cd ${AIE_BUILD_DIR}; \
	v++ -c --mode aie --platform=${PLATFORM} -t ${AIE_TARGET} \
		--config ${PROJECT_DIR}/design/aie/aiecompiler.cfg \
		--work_dir=${PROJECT_DIR}/${BUILD_DIR}/Work \
		--include="${PROJECT_DIR}/design/aie" \
		--include=$(DSPLIB_VITIS)/dsp/L1/src/aie \
		--include=$(DSPLIB_VITIS)/dsp/L1/include/aie \
		--include=$(DSPLIB_VITIS)/dsp/L2/include/aie \
		${PROJECT_DIR}/design/aie/graph.cpp
	@echo ""
	@echo "AIE Application Built..."
	@echo "####################################"
	@echo ""

# Building the host executable
${HOST_BUILD_DIR}/${HOST_EXE}: ${AIE_BUILD_DIR}/libadf.a ${BUILD_DIR}/Work/ps/c_rts/aie_control.cpp design/host/* ${PROJECT_DIR}/design/common.h
	mkdir -p ${HOST_BUILD_DIR}; \
	cd ${HOST_BUILD_DIR}; \
	$(MAKE) -C ${PROJECT_DIR}/design/host/ BUILD_DIR=${PROJECT_DIR}/${BUILD_DIR} -B
	@echo ""
	@echo "Host Application Built..."
	@echo "####################################"
	@echo ""

# Building the XSA
${XSA_BUILD_DIR}/${XSA}: ${PL_BUILD_DIR}/dma_pkt_router.xo ${AIE_BUILD_DIR}/libadf.a design/system_cfgs/system.cfg
	mkdir -p ${XSA_BUILD_DIR}; \
    cd ${XSA_BUILD_DIR}; \
	v++ -g -l --platform ${PLATFORM} -t ${TARGET} \
	    --save-temps \
		--verbose \
		--config ${PROJECT_DIR}/design/system_cfgs/system.cfg \
		-o ${XSA} \
		${PROJECT_DIR}/${AIE_BUILD_DIR}/libadf.a \
		${PROJECT_DIR}/${PL_BUILD_DIR}/dma_pkt_router.xo
	@echo ""
	@echo "XSA Built..."
	@echo "####################################"
	@echo ""

# Running aiesimulator to get print statements and aie.vcd file
${AIESIM_BUILD_DIR}/aie.vcd: ${AIE_BUILD_DIR}/libadf.a
	mkdir -p ${AIESIM_BUILD_DIR}; \
	cd ${AIESIM_BUILD_DIR}; \
	aiesimulator --profile --pkg-dir=${PROJECT_DIR}/${BUILD_DIR}/Work \
		--dump-vcd aie \
		--input-dir ${PROJECT_DIR}/${PLSIM_BUILD_DIR}/plsimulator_output 2>&1 | tee aiesim.log
	@echo ""
	@echo "AIE Simulation, With Profiling, Complete..."
	@echo "####################################"
	@echo ""


# Clean target
clean:
	rm -rf ${BUILD_DIR}
	$(MAKE) -C design/host clean

