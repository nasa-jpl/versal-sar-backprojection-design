# By: Austin Owens
# Date: 9/3/2024
# Desc: Helper script to create custom dts for Versal from xsa
#
# You can run this tcl script by doing the following in 
# versal-design-build/helper_scripts:
# source env_setup.sh && xsct create_dtc.tcl
#
# After running the tcl script, you will see the dts located here:
# custom-device-tree/custom-device-tree/vck190/psv_cortexa72_0/device_tree_domain/bsp/system.dts
#
# In order to register the dts in the yocto build, rename and place 
# the system.dts file within the yocto project under 
# `versal-yocto-build/workspace/poky/meta-jpl-versal/recipes-bsp/device-tree/files/versal.dts`
# and re-run the yocto build and re-create the SDK (as shown in the 
# `versal-manifest` README.md).

set design_build_path [file dirname [file dirname [file normalize [info script]]]]

createdts -zocl -hw $design_build_path/build/hw/xsa/_x/link/int/vpl_gen_fixed.xsa -platform-name vck190 -git-branch xlnx_rel_v2024.1 -compile -out ./custom-device-tree
