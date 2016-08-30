# -*- coding: utf-8 -*-
"""
Created on Thu Jul 28 16:32:47 2016

@author: anindya

Recovery and Relevance scores from comparing Biclustering results against true Biclusters. Works with bibench.

"""

import bibench.all as bb

import sys

def print_score(input_data,input_expected,resultfile):
    expected=bb.read_biclusters(input_expected)
    bcs=bb.read_biclusters(input_data)
    score=bb.jaccard_list(expected,bcs)
    outfile=open(resultfile,"a")
    outfile.write("%s\t%0.2f\t%0.2f\n" %(input_data,score[0],score[1]))
    outfile.close()


input_data=sys.argv[1] # Biclustering results
input_expected=sys.argv[2] # True Biclusters
resultfile=sys.argv[3] # Output scores

print_score(input_data,input_expected,resultfile)


