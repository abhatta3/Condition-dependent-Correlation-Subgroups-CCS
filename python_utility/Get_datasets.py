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
data_const2, expected_const2 = bb.make_const_data(nrows=500, ncols=100, nclusts=1, bicluster_signals=5, shuffle=True)

data_shift, expected_shift = bb.make_shift_data(nrows=100, ncols=100, nclusts=1, shuffle=True)
data_scale, expected_scale = bb.make_scale_data(nrows=100, ncols=100, nclusts=1, shuffle=True)
data_shift_scale, expected_shift_scale = bb.make_shift_scale_data(nrows=100, ncols=100, nclusts=1, shuffle=True)



print_data(bb,data_const1, expected_const1,"Constant_A")
print_data(bb,data_const2, expected_const2,"Constant_B")


print_data(bb,data_shift,expected_shift,"Shift")
print_data(bb,data_scale,expected_scale,"Scale")
print_data(bb,data_shift_scale,expected_shift_scale,"Shift_Scale")

 

