# Versal SAR Backprojection Design

## Introduction

**NOTE: This repository should not be cloned directly! Instead, the following
manifest repo should be used to collect all the necessary repos required for 
building this design: https://github.jpl.nasa.gov/atowens/versal-manifest. Read 
the README.md in that repo first if you want to build this design for the Versal.**

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

The various helper scripts in `versal-design-build` repo are summarized below: 

- `1_copy_tftpboot_files.sh` - Copies TFTP boot files from the Yocto build into your tftpboot directory. Adjust script as needed for your environment.
- `2_copy_nfs_files` - Untars the rootfs that was generated from the Yocto build into your NFS directory. Adjust script as needed for your environment.
- `3_copy_app_files` - Copies the application design files into your NFS directory at `/home/root/app`. Adjust script as needed for your environment.
- `create_dts.tcl` - Create a starting DTS template for Versal built from your generated XSA in the Versal application design project. Adjust DTS as needed.
- `env_setup.sh` - Script to source to set enviornment variables for your system. Adjust script as needed for your enviornment.
- `flash_bootbin_xsct.tcl` - Flashes the `BOOT.BIN` generated from the Versal application design project over JTAG to Versal.

## Flashing Project to Versal Hardware

Although it is possible to place everything onto an SD card, the recommended way to put the OS and application on the hardware during development is over JTAG
and Ethernet using Tiny File Transfer Protocol (TFTP) and Network File System (NFS). This gives the user the most amount of control when configuring the system.

### Tiny File Transfer Protocol (TFTP)

### Network File System (NFS)

NFS facilitates transfer of file over the network. In this case, we are using it to transfer the entire rootfs created by the Yocto build over the network to the 
Linux kernel running on the Versal. If you do not have NFS on your host machine, you can install it with `sudo apt install nfs-kernel-server`.

You will need to modify your `/etc/exports` file. Depending on your desired NFS configuration and IP address, you may have a diffrent setup. An example of the setup
could look something like this:

```console
/srv/nfs 192.168.2.0/24(rw,no_root_squash,no_all_squash,crossmnt)
```

After you edit the `/etc/exports` file, you need to make the file share available and restart the NFS with the following commands:

```bash
sudo exportfs -a
sudo systemctl restart nfs-kernel-server
```

**WARNING:** The Linux kernel running on the Versal must have its kernel arguments point to the NFS share. There are multiple ways to update the Linux kernel's bootargs. The 
recommended way is through pxebooting via TFTP. Refer to the TFTP section above.

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
ACAP are stored in their own repositories. For example, the build and helper
scripts are often repetitive between Versal designs, so they are pulled out
into their own repository 
([versal-design-build](https://github.jpl.nasa.gov/atowens/versal-design-build)).
All the repos necessary for building the Linux OS running on the ARM of the
Versal is also pulled out into their own repository (see the manifest repo
linked at the top of this README.md file for more info). This allows for
easier management and upgrades to various elements of Versal projects since only
one repository needs to be updated for all Versal designs to get those features.

Of course, like any organizational structure, there will also be negatives.
For example, updates made to a shared repo could potentially break other Versal
designs. To combat this, there currently are AIE sims that can be ran to verify
parts of the design. The level of regression testing is currently in it's
infancy and should be further matured to be more comprehensive and automated,
such as having sims for the host code that runs within the ARM on the Versal,
or creating a CI-CD pipeline to automatically run regression tests with AMD's
emulation tools when any changes are made to the repo.
