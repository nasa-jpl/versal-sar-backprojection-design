#!/usr/bin/python3
# By: Austin Owens
# Date: 6/3/2024
# Desc: Creates FFT input data in 128 bit chunks for AIE


import os

mat_dims = [16, 32, 64, 128, 256, 1024]
iter_cnt = 1

os.system("rm -rf c*")

for mat_dim in mat_dims:
    mat_dir = "cint16_simple_mat"
    os.system("mkdir -p " + mat_dir)

    cint16_mat_file = open("{0}/mat_{1}x{1}.txt".format(mat_dir, mat_dim), "a+")
    for j in range(0, iter_cnt):
        for i in range(0, mat_dim * mat_dim):
            if i == 0:
                cint16_mat_file.write("1 1" + "\n")
            else:
                cint16_mat_file.write("0 0" + "\n")
cint16_mat_file.close()

for mat_dim in mat_dims:
    mat_dir = "cfloat_simple_mat"
    os.system("mkdir -p " + mat_dir)

    cfloat_mat_file = open("{0}/mat_{1}x{1}.txt".format(mat_dir, mat_dim), "a+")
    for j in range(0, iter_cnt):
        for i in range(0, mat_dim * mat_dim):
            if i == 0:
                cfloat_mat_file.write("1.5 1.5" + "\n")
            else:
                cfloat_mat_file.write("0 0" + "\n")
cfloat_mat_file.close()

os.system("chmod 755 -R *")
