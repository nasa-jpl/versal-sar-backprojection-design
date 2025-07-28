#!/bin/bash
#
# By: Austin Owens
# Date: 5/21/2024
# Desc: Source this script to allow the makefiles to build the system properly

WORKSPACE_PATH=$(dirname $(dirname $(dirname $(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd))))

# =======================================================
# Set Platform ,Vitis Versal Image, SDK, and PDM repos
# =======================================================
export PLATFORM_REPO_PATHS="/mnt/disk/Xilinx/Vitis/2024.1/base_platforms"
export XILINX_VITIS="/mnt/disk/Xilinx/Vitis/2024.1"
export PDM_PATH="/mnt/disk/Xilinx/PDM/2024.1.1"
export SDK_PATH="$WORKSPACE_PATH/versal-yocto-build/workspace/poky/build/tmp/deploy/sdk/toolchain"

# =======================================================
# Set ARM DTB, bootloader, Linux kernel, and rootfs paths
# =======================================================
export DTB="$WORKSPACE_PATH/versal-yocto-build/workspace/poky/build/tmp/deploy/images/vck190-versal/system.dtb"
export BL31_ELF="$WORKSPACE_PATH/versal-yocto-build/workspace/poky/build/tmp/deploy/images/vck190-versal/arm-trusted-firmware.elf"
export UBOOT="$WORKSPACE_PATH/versal-yocto-build/workspace/poky/build/tmp/deploy/images/vck190-versal/u-boot.elf"
export IMAGE="$WORKSPACE_PATH/versal-yocto-build/workspace/poky/build/tmp/deploy/images/vck190-versal/Image"
export ROOTFS="$WORKSPACE_PATH/versal-yocto-build/workspace/poky/build/tmp/deploy/images/vck190-versal/jpl-versal-image-vck190-versal.rootfs.tar.gz"

# ========================================================
# Set DSP Library for Vitis
# ========================================================
export DSPLIB_VITIS="/mnt/disk/Xilinx/Vitis_Libraries"

# ====================================================
# Source Versal Image ,Vitis and Aietools
# ====================================================
# Run the below command to setup environment and CXX
source $SDK_PATH/environment-setup-cortexa72-cortexa53-poky-linux
source $XILINX_VITIS/settings64.sh

# ====================================================
# Source Versal Image ,Vitis and Aietools
# ====================================================
# Run the below command to setup environment for power estimation
source $PDM_PATH/settings64.sh

# =========================================================
# Platform Selection...
# =========================================================
tgt_plat=xilinx_vck190_base_202410_1
#tgt_plat=base_pfm_adm_pa100
export PLATFORM=$PLATFORM_REPO_PATHS/$tgt_plat/$tgt_plat\.xpfm

# ==========================================================
# Validating Tool Installation
# ==========================================================
echo ""
echo "Aiecompiler:"
which aiecompiler
echo ""
echo "Vivado:"
which vivado
echo ""
echo "Vitis:"
which vitis
echo ""
echo "Vitis HLS:"
which vitis_hls
echo ""
echo "DSPLIBS"
echo "$DSPLIB_VITIS"
echo ""
echo "SDK"
echo "$SDK_PATH"
echo ""
echo "PDM"
echo "$PDM_PATH"
echo ""
