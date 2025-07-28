#!/bin/bash
#
# By: Austin Owens
# Date: 9/3/2024
# Desc: Copies the application design files into the /home/root/app dir.
# Adjust as needed for your environment.

WORKSPACE_PATH=$(dirname $(dirname $(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)))
PACKAGE_PATH=${WORKSPACE_PATH}/versal-design-build/build/hw/package
NFS_PATH=/nfs/versal/rootfs

echo "BEFORE (/home/root/app):"
sudo ls -l ${NFS_PATH}/home/root/app 2>/dev/null

# Remove app dir
sudo rm -rf ${NFS_PATH}/home/root/app

# Make app directory
sudo mkdir -p ${NFS_PATH}/home/root/app

# Copy over a.xclbin, application, and run script
sudo cp -r ${PACKAGE_PATH}/sd_card/a.xclbin ${PACKAGE_PATH}/sd_card/sar_backproject.elf ${PACKAGE_PATH}/sd_card/run_script_hw.sh ${NFS_PATH}/home/root/app/

# Copy over slowtime dataset file
sudo cp -r ${PACKAGE_PATH}/sd_card/gotcha_slowtime_pass1_360deg_HH.csv ${NFS_PATH}/home/root/app/

# Copy over range compression dataset file
RC_SAMPLES=$(grep '^#define RC_SAMPLES' ${WORKSPACE_PATH}/versal-design-build/design/common.h | awk '{print $3}')
sudo cp -r ${PACKAGE_PATH}/sd_card/gotcha_${RC_SAMPLES}-out-of-424-rc-samples_pass1_360deg_HH.csv ${NFS_PATH}/home/root/app/

# Make permissions on application and runscript executable
sudo chmod +x ${NFS_PATH}/home/root/app/sar_backproject.elf ${NFS_PATH}/home/root/app/run_script_hw.sh

echo "AFTER (/home/root/app):"
sudo ls -l ${NFS_PATH}/home/root/app
