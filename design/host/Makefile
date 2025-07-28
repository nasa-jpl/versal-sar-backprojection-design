# By: Austin Owens
# Date: 5/21/2024
# Desc: Makefile for building the host application that runs on the ARM processor

SDKTARGETSYSROOT ?= ${SYSROOT}
OUTPUT_DIR = ${BUILD_DIR}/host

HOST_EXE = sar_backproject.elf
HOST_OBJ = ${OUTPUT_DIR}/aie_control_xrt.o ${OUTPUT_DIR}/sar_backproject.o ${OUTPUT_DIR}/main.o 
AIE_CTRL_CPP = ${BUILD_DIR}/Work/ps/c_rts/aie_control_xrt.cpp


# __AIE_ARCH__ is specified in <Vitis_DSP_Lib>/Vitis_Libraries/dsp/L1/include/aie/device_defs.h
# The following values are 10 (AIE1), 20 (AIE2, aka AIE-ML), 21 (AIE_PS), 22 (AIE_PS2)
CXXFLAGS += -std=c++17 -O0 -g -Wall -c -fmessage-length=0 --sysroot=${SDKTARGETSYSROOT} \
			-D__AIE_ARCH__=10 \
		    -fopenmp \
			-I$(XILINX_VIVADO)/include/ \
			-I${SDKTARGETSYSROOT}/usr/include/xrt/ \
			-I${XILINX_VITIS}/aietools/include \
			-I$(DSPLIB_VITIS)/dsp/L1/include/aie \
			-I./ \
			${HOST_INC}

LDFLAGS += --sysroot=${SDKTARGETSYSROOT} \
		   -ladf_api_xrt \
		   -lgcc \
		   -lc \
		   -lxrt_coreutil \
		   -lxilinxopencl \
		   -lpthread \
		   -lrt \
		   -ldl \
		   -lcrypt \
		   -lstdc++ \
		   -lhdf5 \
		   -fopenmp \
		   -L${SDKTARGETSYSROOT}/usr/lib/ \
		   -L$(XILINX_VITIS)/aietools/lib/aarch64.o

# Use OUTPUT_DIR for the output directory
${OUTPUT_DIR}/${HOST_EXE}: ${HOST_OBJ}
	${CXX} -o $@ $^ ${LDFLAGS}

${OUTPUT_DIR}/%.o: %.cpp
	mkdir -p ${OUTPUT_DIR}
	${CXX} ${CXXFLAGS} -o $@ $<

${OUTPUT_DIR}/main.o: ./*  # Rebuild main.o if anything in SRC_DIR changes
	mkdir -p ${OUTPUT_DIR}
	${CXX} ${CXXFLAGS} -c -o $@ ./main.cpp

${OUTPUT_DIR}/aie_control_xrt.cpp: ${AIE_CTRL_CPP}
	cp -f ${AIE_CTRL_CPP} ${OUTPUT_DIR}/aie_control_xrt.cpp

clean:
	rm -rf ${OUTPUT_DIR}/*.o .Xil ${OUTPUT_DIR}/${HOST_EXE} ${OUTPUT_DIR}/aie_control_xrt.cpp

