# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:32:47 2016

@author: anindya

Datasets for Biclustering

"""

import bibench.all as bb


def print_data(bb,data,expected,string):
    bfile="Data_"+string+"_bicluster.txt"
    bb.write_expression_data(data, bfile, sep='\t', genes=None, conditions=None)
    tfile="True_"+string+"_bicluster.txt"
    bb.write_biclusters(expected,tfile)



data_const1, expected_const1 = bb.make_const_data(nrows=100, ncols=100, nclusts=1, bicluster_signals=5, shuffle=True)
data_const2, expected_const2 = bb.make_const_data(nrows=250, ncols=100, nclusts=1, bicluster_signals=5, shuffle=True)

data_shift_scale1, expected_shift_scale1 = bb.make_shift_scale_data(nrows=100, ncols=100, nclusts=1, shuffle=True)
data_shift_scale2, expected_shift_scale2 = bb.make_shift_scale_data(nrows=250, ncols=100, nclusts=1, shuffle=True)
data_shift_scale4, expected_shift_scale4 = bb.make_shift_scale_data(nrows=250, ncols=100, nclusts=3, shuffle=True)


print_data(bb,data_const1, expected_const1,"Constant_100_1")
print_data(bb,data_const2, expected_const2,"Constant_250_1")


print_data(bb,data_shift_scale1,expected_shift_scale1,"Shift_Scale_100_1")
print_data(bb,data_shift_scale2,expected_shift_scale2,"Shift_Scale_250_1")
print_data(bb,data_shift_scale4,expected_shift_scale4,"Shift_Scale_250_3")

 

