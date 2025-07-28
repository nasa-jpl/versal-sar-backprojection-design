#!/bin/bash
#
# By: Austin Owens
# Date: 9/3/2024
# Desc: Untars the rootfs that was generated from the Yocto build into the NFS dir.
# Adjust as needed for your environment.

WORKSPACE_PATH=$(dirname $(dirname $(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &>/dev/null && pwd)))
NFS_PATH=/nfs/versal/rootfs

echo "BEFORE:"
ls -l ${NFS_PATH}

sudo rm -r ${NFS_PATH}/*
sudo tar -xvf ${WORKSPACE_PATH}/versal-yocto-build/workspace/poky/build/tmp/deploy/images/vck190-versal/jpl-versal-image-vck190-versal.rootfs.tar.gz -C ${NFS_PATH}

echo "AFTER:"
ls -l ${NFS_PATH}
