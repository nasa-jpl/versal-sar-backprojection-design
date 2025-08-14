# Versal SAR Backprojection Design

## Introduction

**NOTE: This repository should not be cloned directly! Instead, follow the instructions in the 
[versal-manifest](https://github.jpl.nasa.gov/atowens/versal-manifest) repo to collect all the 
necessary repos required for building this design.**

This repo contains the source code for implementing backprojection onto 
the AMD Versal ACAP, along with a makefile to help build the project and 
various helper scripts to help flash the project for the Versal target.

The repo folder hierarchy contains the following:

Directory                     | Contents
------------------------------|-----------
design/aie                    | Source code for AIE engine on Versal
design/exec_scripts           | Scripts that gets transferred to Versal ARM to execute the design
design/host                   | Source code for the ARM Cortex-A72 on Versal
design/pl                     | Source code for the FPGA on Versal (currently not being used)
design/profiling_cfgs         | Config that gets transferred to Versal ARM describing various setting of how the design should be ran
design/system_cfgs            | Provides information to compiler when compiling the entire design
design/test_data              | Test data for slowtime and range compression samples
design/vivado_metrics_scripts | Scripts for executing various metrics (power, resource utilization, etc.) for compiled design
design/common.h               | Configurations for the design used throughout the entire application
helper_scripts                | Helper scripts used to help setup your environment variables and flash the applications on the Versal
Makefile                      | Makefile to build out entire project
doc                           | Design documentation

### Git LFS Installation

This repo uses Git LFS for it's data input files. In order to use Git LFS, you'll
need to download and install a new program that's separate from Git.

1. Install git-lfs via your package manager. Ex: `sudo apt install git-lfs`. 
Alternatively, you can navigate to [git-lfs.com](https://git-lfs.com/) and 
click Download.
2. `cd versal-sar-backprojection-design`
3. `git lfs install`
4. `git lfs pull`

### Helper Scripts

The various helper scripts in `versal-sar-backprojection-design` repo are summarized below: 

- `1_copy_tftpboot_files.sh` - Copies TFTP boot files from the Yocto build into your tftpboot directory. Adjust script as needed for your environment.
- `2_copy_nfs_files` - Untars the rootfs that was generated from the Yocto build into your NFS directory. Adjust script as needed for your environment.
- `3_copy_app_files` - Copies the application design files into your NFS directory at `/home/root/app`. Adjust script as needed for your environment.
- `create_dts.tcl` - Create a starting DTS template for Versal built from your generated XSA in the Versal application design project. Adjust DTS as needed.
- `env_setup.sh` - Script to source to set environment variables for your system. Adjust script as needed for your environment.
- `flash_bootbin_xsct.tcl` - Flashes the `BOOT.BIN` generated from the Versal application design project over JTAG to Versal.

## Building the Project

After fetching the required dependencies for this design (as shown in the 
[versal-manifest](https://github.jpl.nasa.gov/atowens/versal-manifest) repo), the different components of the design 
can be built with the below instructions.

### Yocto Build

From your workspace root directory:

```bash
cd versal-sar-backprojection-design
docker build --build-arg USERNAME=<user_name> --build-arg UID="$(id -u)" -t yocto-img .
./docker.run -u <user_name> -s $(pwd)/workspace -i yocto-img /bin/bash

# Use <user_name> for the password below
sudo chown -R <user_name> workspace

source build.sh
```

Where `<user_name>` is any user name you want to have for inside the container. This 
does not need to be your host computer's username.

The Yocto build can take an hour or longer to build.

### SAR Application Design Build

If you just built the Yocto build from the above step, exit out of the docker container
with `exit`. 

You must export the Yocto SDK so that the design can utilize the correct toolchain and 
target OS's sysroot. To do this, execute the following from your workspace root directory:

```bash
cd versal-yocto-build/workspace/poky/build/tmp/deploy/sdk
./poky-glibc-x86_64-jpl-versal-image-cortexa72-cortexa53-vck190-versal-toolchain-5.0.3.sh -d toolchain -y
```

The above command install Yocto's SDK into a directory called `toolchain` in the current directory. 
This allows the design application code to link against the sysroot that Yocto created.

Next, you need to make sure the [helper/scripts/env_setup.sh](helper_scripts/env_setup.sh) contains the 
correct paths to your Xilinx installation.

Lastly, build the design application:

```bash
cd versal-sar-backprojection-design
source helper_scripts/env_setup.sh
make
```

NOTE: You can run `make` with additional parameters. Refer to top comment in the [Makefile](Makefile) to
see what options are available.

## Flashing Project to Versal Hardware

At this point, you should have everything built from the previous section.

Although it is possible to place everything onto an SD card, the recommended way to put the OS and application on the hardware during development is over JTAG
and Ethernet using Tiny File Transfer Protocol (TFTP) and Network File System (NFS). This gives the user the most amount of control when configuring the system.

### Network File System (NFS)

NFS facilitates transfer of file over the network. In this case, we are using it to transfer the entire rootfs created by the Yocto build over the network to the 
Linux kernel running on the Versal. If you do not have NFS on your host machine, you can install it with `sudo apt install nfs-kernel-server`.

You will need to modify your `/etc/exports` file. Depending on your desired NFS configuration and IP address, you may have a diffrent setup. An example of the setup
could look something like this:

```console
/nfs/versal/rootfs 192.168.10.0/24(rw,no_root_squash,no_all_squash,crossmnt)
```

After you edit the `/etc/exports` file, you need to make the file share available and restart the NFS with the following commands:

```bash
sudo exportfs -a
sudo systemctl restart nfs-kernel-server
```

The [helper_scripts/2_copy_nfs_files.sh](helper_scripts/2_copy_nfs_files.sh) untars the rootfs that was generated from the Yocto build into the NFS directory.
**YOU MUST ADJUST THIS FILE AS NEEDED FOR YOUR ENVIRONMENT.** The `2_copy_nfs_files.sh` assumes your NFS directory is `/nfs/versal/rootfs` by default.
If you have followed the above NFS configuration, you shouldn't have to modify this file. However, if you have a diffrent NFS configuration, change the
`2_copy_nfs_files.sh` script accordingly.

**IMPORTANT:** The Linux kernel running on the Versal must have its kernel arguments (bootargs) point to the NFS share running on your host computer. 
There are multiple ways to update the Linux kernel's bootargs. The recommended way is through the pxeboot configuration located at 
[helper_scripts/pxelinux.cfg/default-arm-versal](helper_scripts/pxelinux.cfg/default-arm-versal). By default, the following boot arguments for the `nfsroot` and 
`ip` args are:
- `nfsroot=192.168.10.1:/nfs/versal/rootfs,tcp,nfsvers=3`
    - `192.168.10.1` - NFS server IP (your host computer)
    - `/nfs/versal/rootfs` - Path to the NFS directory containing the rootfs.
    - `tcp` - Use TCP transport for NFS (more reliable than UDP).
    - `nfsvers=3` - Force NFSv3 (common for embedded roots)
- `ip=192.168.10.2:::255.255.255.0:versal:eth0:bootp`
    - `192.168.10.2` - Versal IP
    - `255.255.255.0` - Netmask
    - `versal` - Hostname the kernel will assign
    - `eth0` - Network interface to bring up early (for fetching NFS root).
    - `bootp` - Autoconfiguration method (accepts off|on|dhcp|bootp|rarp; here it's BOOTP/DHCP-style).

Change the above arguments as appropriate within the `default-arm-versal` file for your networking configuration. The TFTP client (covered in the next section)
is responsible for transfering over the `default-arm-versal` file to U-Boot so it can utilize these bootargs when launcing the Linux kernel.

### Tiny File Transfer Protocol (TFTP)

TFTP is responsible for transferring the Linux Kernel and [pxeboot file (default-arm-versal)](helper_scripts/pxelinux.cfg/default-arm-versal) from
your host machine over to the Versal. This README will not cover TFTP installation and configuration since there are several online tutorials that cover this process.

The [helper_scripts/1_copy_tftpboot_files.sh](helper_scripts/1_copy_tftpboot_files.sh) copies the TFTP boot files from the Yocto build into your host computer's 
tftpboot directory. **YOU MUST ADJUST THIS FILE AS NEEDED FOR YOUR ENVIRONMENT.** The `1_copy_tftpboot_files.sh` assumes your tftp boot directory is `/tftpboot`
by default. If this is not the case, change the  `1_copy_tftpboot_files.sh` script accordingly.

### Application Executable

Once the NFS is setup and you have the design built, you can use the [helper_script/3_copy_app_files.sh](helper_script/3_copy_app_files.s) script to copy the application
design files into the `/home/root/app` directory on the Versal/host NFS. **YOU MUST ADJUST THIS FILE AS NEEDED FOR YOUR ENVIRONMENT.** The `3_copy_app_files.sh` assumes
your NFS directory is `/nfs/versal/rootfs` by default. If you have followed the above NFS configuration, you shouldn't have to modify this file. However, if you have a
diffrent NFS configuration, change the `3_copy_app_files.sh` script accordingly.

### Serial Console

Before flashing the design over JTAG, open up a serial console session with the Versal. The serial console should be accessible via `/dev/ttyUSB3` (although the
serial devices are dynamically assigned and can change). Both minicom and screen can be used to access it. 

For example:

```console
minicom -D /dev/ttyUSB3

# OR

screen /dev/ttyUSB3 115200
# To quit screen, press `ctrl+a` then `k` to quit. To detach from screen press `ctrl+a` then `d`. To re-attach to running screen run: `screen -r`
```

The "System Controller" also has a serial console on a separate port that will start printing boot messages immediately upon power-up. When finished booting, it will print
a link for accessing the web application called "BEAM".

Upon the first time flashing over JTAG, all serial console ports should be displayed in order to know which serial device port belongs to the Linux application that will boot.

### JTAG

The [helper_scripts/flash_bootbin_xsct.tcl](helper_scripts/flash_bootbin_xsct.tcl) script can be used to flash `BOOT.BIN` into Versal over JTAG. The following command to execute
this TCL script is the following:

```console
cd versal-sar-backprojection-design
xsct ./helper_scripts/flash_bootbin_xsct.tcl
```

Check the serial console to see JTAG/U-Boot/Linux messages.

IMPORTANT: If you are not using a DHCP server for within your network configuration, you will need to interrupt U-Boot in the serial console by pressing "enter" (or `ctrl+c`) while
it is trying to boot. Once you have a U-Boot shell, enter in the following network settings:

```console
setenv ipaddr <versal_ip>
setenv netmask <netmask>
setenv serverip <nfs_server_ip>
ping <nfs_server_ip> # OPTIONAL: Make sure you can ping host
setenv bootcmd 'if pxe get; then pxe boot; fi'
boot
```
Where:
- `<versal_ip>`: IP address of the Versal
- `<netmask>`: Your netmask IP configuration (e.g. `255.255.255.0`)
- `<nfs_server_ip>`: IP address of host NFS server

The above network settings should match the network settings within your customized [pxeboot configuration file](helper_scripts/pxelinux.cfg/default-arm-versal).

## Documentation Generation

The following packages were installed for documentation generation:
- `latexmk` - Automates compiling LaTeX documents
- `zathura` - Document viewer
- `zathura-pdf-poppler` - Plugin to add PDF support for Zathura
- `texlive-latex-extra` - Adds .svg image support
- `inkscape` - Graphics support

## FAQ

1. Why doesn't this repo contain everything necessary to build this design?

In order to promote better version control practices, dependency management,
and extensibility, different elements of building out designs for the Versal
ACAP are stored in their own repositories. For example, All the repos necessary
for building the Linux OS running on the ARM of the Versal are pulled out into
their own repository (see the [versal-manifest](https://github.jpl.nasa.gov/atowens/versal-manifest)
repo to see other repo dependancies needed to build/flash this design to a Versal).
This allows for easier management and upgrades to various elements of Versal projects
since only one repository needs to be updated for all Versal designs to get those features.

Like any organizational structure, there will also be negatives. For example,
updates made to a shared repo could potentially break other Versal designs.
To combat this, there currently are AIE sims that can be ran to verify parts
of the design. The level of regression testing is currently in it's infancy
and should be further matured to be more comprehensive and automated, such
as creating a CI-CD pipeline to automatically run regression tests with AMD's
emulation tools when any changes are made to the repo.
