# Versal SAR Backprojection Design

## Introduction

**NOTE: This repository should not be cloned directly! Instead, the following
manifest repo should be used to collect all the necessary repos required for 
building this design: https://github.jpl.nasa.gov/atowens/versal-manifest. Read 
the README.md in that repo for further instructions.**

This repo only contains the source code for implementing backprojection onto 
the AMD Versal ACAP and does not contain various helper scripts or makefiles
to help build the project for the Versal target.

## Git LFS Installation

This repo uses Git LFS for it's documentation. In order to use Git LFS, you'll
need to download and install a new program that's separate from Git.

1. Install git-lfs via your package manager. Ex: `sudo apt install git-lfs`. 
Alternatively, you can navigate to [git-lfs.com](https://git-lfs.com/) and 
click Download.
2. `cd versal-design-build/design`
3. `git lfs install`
