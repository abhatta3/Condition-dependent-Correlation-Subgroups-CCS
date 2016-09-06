#!/usr/bin/env python
"""
About: Draw three Bicluster graphs for testing condition dependency. 

Dependency: matplotlib, networkx 

Execution:  python Condition_dependency_graph.py  argument_1  argument_2

argument_1:  is the input data matrix used for biclustering.

argument_2:  is an input file with two biclusters.


Example: python Condition_dependency_graph.py GDS589_full.soft Bicluster_network_input.txt 




Output:

3 output figures:


0.png: Correlation graph plot on all the samples in the input dataset for correlation>0.8 as edges, 

       blue nodes are from bicluster 1, green nodes are from bicluster 2 and white nodes are neighboring nodes with correlation  



1.png: Correlation graph plot on all the bicluster 1 samples for correlation>0.8 as edges, 

       blue nodes are from bicluster 1, green nodes are from bicluster 2 and white nodes are neighboring nodes with correlation  

       

2.png: Correlation graph plot on all the bicluster 2 samples for correlation>0.8 as edges, 

       blue nodes are from bicluster 1, green nodes are from bicluster 2 and white nodes are neighboring nodes with correlation 
       
"""
__author__ = """anb013@eng.ucsd.edu"""


import sys

import networkx as nx
import matplotlib.pyplot as plt

import math

def average(x):
    assert len(x) > 0
    return float(sum(x)) / len(x)

def pearson_def(x, y):
    assert len(x) == len(y)
    n = len(x)
    assert n > 0
    avg_x = average(x)
    avg_y = average(y)
    diffprod = 0
    xdiff2 = 0
    ydiff2 = 0
    for idx in range(n):
        xdiff = x[idx] - avg_x
        ydiff = y[idx] - avg_y
        diffprod += xdiff * ydiff
        xdiff2 += xdiff * xdiff
        ydiff2 += ydiff * ydiff

    return diffprod / math.sqrt(xdiff2 * ydiff2)
    
    


input_data_file=sys.argv[1]

input_bicluster=sys.argv[2]


with open(input_bicluster,"r") as bfile:
    for i,line in enumerate(bfile):
        if i>5:
            break
        if i==1:
            bicluster_1_g=line.strip().split()
        elif i==2:
            bicluster_1_s=line.strip().split()
        elif i==4:
            bicluster_2_g=line.strip().split()
        elif i==5:
            bicluster_2_s=line.strip().split()


genelist=bicluster_1_g+bicluster_2_g 


data_dic={}

with open(input_data_file,"r") as data:
    for i,line in enumerate(data):
        if i==0:
            H=line.strip().split("\t")
            continue
        
        dline=line.strip().split("\t")
        for j,cl in enumerate(dline):
            if j==0:
                row_dic={}
                gene=cl
                data_dic[gene]=row_dic
                continue
            
            sample=H[j]

            data_dic[gene][sample]=float(cl)
            


gene_dic={}
color_dic={}

for item in genelist:
    if item in bicluster_1_g and item in bicluster_2_g:
        color_dic[item]='r'
    elif item in bicluster_1_g:
        color_dic[item]='b'
    elif item in bicluster_2_g:
        color_dic[item]='g'
        
    gene_dic[item]=1



def get_graph(genelist,gene_dic,bicluster_1_g,bicluster_2_g,color_dic,term,samplelist,data_dic,figurname):
   
   
    n_dic={}
    
    edg_dic={}
                
    for g in gene_dic:
        nd={}
        edg_dic[g]=nd
        for ng in data_dic:
            if g==ng:
                continue
            
            if term==1 and g in bicluster_1_g and ng in bicluster_1_g:
                edg_dic[g][ng]='y'
                continue
            elif term==2 and g in bicluster_2_g and ng in bicluster_2_g:
                edg_dic[g][ng]='y'
                continue
            
            gl=[]
            ngl=[]
            
            for s in samplelist:
                gl.append(data_dic[g][s])
                
            for s in samplelist:
                ngl.append(data_dic[ng][s])
            pcr=pearson_def(gl, ngl)    
            if pcr>0.8 or pcr<-0.8:
                edg_dic[g][ng]='y'
                if not gene_dic.has_key(ng):
                    n_dic[ng]=1
                    color_dic[ng]='w'
    
    
    for g in n_dic:
        gene_dic[g]=1
        nd={}
        edg_dic[g]=nd
        
        
    
    edl=[]
    edcolor=[]
    node_c=[]
    node_l=[]
    labels={}
    count=0
    
    
    selfn={}
    
    for i,g in enumerate(gene_dic):
        selfn[g]=1
        for j,ng in enumerate(gene_dic):
            if j<=i:
                continue
            if edg_dic[g].has_key(ng):
                edl.append((i,j))
                selfn[g]=0
                selfn[g]=0
                edcolor.append(edg_dic[g][ng])
                if i not in node_l:
                    node_l.append(i)
                    ncv=color_dic[g]
                    node_c.append(ncv)
                    labels[count]=g
                    count+=1
                if j not in node_l:
                    node_l.append(j)
                    ncv=color_dic[ng]
                    node_c.append(ncv)
                    labels[count]=ng
                    count+=1
    
    for i,g in enumerate(gene_dic):
        if selfn[g]==1 and not n_dic.has_key(g):
            edl.append((i,i))
            edcolor.append('b')
            node_l.append(i)
            ncv=color_dic[g]
            node_c.append(ncv)
            labels[count]=g
            count+=1
    
    G=nx.Graph()
    
    for edge in edl:
        G.add_edge(edge[0], edge[1])
    #pos=nx.shell_layout(G)
    pos=nx.spring_layout(G)
    
    """
    if graph_layout == 'spring':
        pos=nx.spring_layout(G)
    elif graph_layout == 'spectral':
        pos=nx.spectral_layout(G)
    elif graph_layout == 'random':
        pos=nx.random_layout(G)
    else:
        pos=nx.shell_layout(G)
    """
    
    #pos=nx.spring_layout(G) # positions for all nodes
    
    # nodes
    nx.draw_networkx_nodes(G,pos,
                           nodelist=node_l,
                           node_color=node_c,
                           node_size=50,
                       alpha=0.8)
    
    # edges
    #nx.draw_networkx_edges(G,pos,width=1.0,alpha=0.3)
    nx.draw_networkx_edges(G,pos,
                           edgelist=edl,
                           width=1,alpha=0.3,edge_color=edcolor)
    
    
    #nx.draw_networkx_labels(G,pos,labels,font_size=8)
    
    plt.axis('off')
    
    plt.savefig(figurname) # save as png





get_graph(genelist,gene_dic,bicluster_1_g,bicluster_2_g,color_dic,0,H[1:],data_dic,"00.png")
get_graph(genelist,gene_dic,bicluster_1_g,bicluster_2_g,color_dic,1,bicluster_1_s,data_dic,"01.png")
get_graph(genelist,gene_dic,bicluster_1_g,bicluster_2_g,color_dic,2,bicluster_2_s,data_dic,"02.png") 
