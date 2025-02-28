# Versal SAR Backprojection Design

## Introduction

**NOTE: This repository should not be cloned directly! Instead, the following
manifest repo should be used to collect all the necessary repos required for 
building this design: https://github.jpl.nasa.gov/atowens/versal-manifest. Read 
the README.md in that repo for further instructions.**

This repo only contains the source code for implementing backprojection onto 
the AMD Versal ACAP and does not contain various helper scripts or makefiles
to help build the project for the Versal target. 

The repo folder hierarchy contains the following:

Directory              | Contents
-----------------------|-----------
aie                    | Source code for AIE engine on Versal
doc                    | Design documentation
exec_scripts           | Scripts that gets transferred to Versal ARM to execute the design
host                   | Source code for the ARM Cortex-A72 on Versal
pl                     | Source code for the FPGA on Versal (currently not being used)
profiling_cfgs         | Config that gets transferred to Versal ARM describing various setting of how the design should be ran
system_cfgs            | Provides information to compiler when compiling the entire design
vivado_metrics_scripts | Scripts for executing various metrics (power, resource utilization, etc.) for compiled design


## Git LFS Installation

This repo uses Git LFS for it's documentation. In order to use Git LFS, you'll
need to download and install a new program that's separate from Git.

1. Install git-lfs via your package manager. Ex: `sudo apt install git-lfs`. 
Alternatively, you can navigate to [git-lfs.com](https://git-lfs.com/) and 
click Download.
2. `cd versal-design-build/design`
3. `git lfs install`

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
