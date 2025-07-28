#!/bin/bash
#
# By: Austin Owens
# Date: 9/3/2024
# Desc: Copies TFTP boot files from the Yocto build into the tftpboot dir.
# Adjust as needed for your environment.

WORKSPACE_PATH=$(dirname $(dirname $(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)))
SRC_PATH=$WORKSPACE_PATH/versal-yocto-build/workspace/poky/build/tmp/deploy/images/vck190-versal
TFTPBOOT_PATH=/tftpboot

echo "BEFORE:"
ls -l ${TFTPBOOT_PATH}

rm -r ${TFTPBOOT_PATH}/*

cp -a ${SRC_PATH}/Image* ${SRC_PATH}/devicetree ${WORKSPACE_PATH}/versal-design-build/helper_scripts/pxelinux.cfg/ ${SRC_PATH}/system.dtb ${TFTPBOOT_PATH}

chmod -R 777 ${TFTPBOOT_PATH}/*

echo "AFTER:"
ls -l ${TFTPBOOT_PATH}
