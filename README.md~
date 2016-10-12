####################################################################################################
###      Condition-dependent Correlation Subgroup (CCS)                                          ###
###----------------------------------------------------------------------------------------------###
###                                                                                              ###
### This file is part of the CCS package for biclustering analysis                               ###
###                                                                                              ###
###----------------------------------------------------------------------------------------------###


=======================
Introduction to CCS
=======================

Condition-dependent Correlation Subgroup (CCS) is a Biclustering algorithm that integrates several important features for developing an effective algorithm 
for comprehensive discovery of functionally coherent biclusters. 

The most important features of CCS are : 

(i) searching and reporting the bicluster in the decreasing order of variability of the base genes, 
(ii) defining biclusters from expression patterns, 
(iii) defining most comprehensive classes of similarity patterns, 
(iv) defining biclusters as condition dependent modules, 
(v) annotating the biclusters with condition-dependency score, 
(vi) merging the overlapping biclusters based on condition dependency and overlap and 
(vii) considering both negative and positive correlation in the same bicluster. 

The algorithm is implemented in C. See the Steps in A and B for compilation and execution of the C code.

The structure of the CCS algorithm is particularly suitable for parallel computing. A CUDA C based GPGPU computing code is included here as a parallel version of the algorithm. 
You need a programmable GPU card and CUDA C complier for compilation and execution of our code. Follow the steps in C and D.

The performance of CCS was compared on synthetic and real gene expression datasets. We also showed that there is an equivalence between CCS biclusters and condition-dependent co-expression network modules. These codes are implemented in Python 2 and work on Python 2.7 or higher. They are made available in "python_utility" directory. 


Synthetic datasets, CCS results, true biclusters, CCS biclustering results on gene expression datasets and the datasets are available in "Results" directory.



============================
Installation and Execution
============================


A. Compilation of C code
====================================

1. Change your current directory to "src".
	$cd src

2. Type "make" to create executable "ccsbc" in the home directory
	$make
3. Back to parent directory
        cd ../

B. Execute C code
====================================

Type following in the linux command line

	./ccsbc -t [correlation threshold] -i [input file] -o [output file]

Parameters:

-t [0-1.0]: Specify correlation threshold "theta" between 0 to 1. Recommended value is 0.8. 
-i [input_data_file]
-o [output_data_file]


Example: $./ccsbc -t 0.8 -i ./Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output.txt


Additional parameters:

-m [1 - number of gene/rows in the data matrix]: Set the number of base gene that are to be considered for forming biclusters. 
 						 Default value is 1000 or maximum number of genes when that is less than 1000. 

Example: 

$./ccsbc -t 0.8 -i ./Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output.txt -m 50



-g [0.0 - 100.0]: Minimum gene set overlap required for merging the overlapped biclusters. 
 		  Default value is 100.0 for 100% overlap. 

Example: 

$./ccsbc -t 0.8 -i ./Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output.txt -g 50.0


-p [0/1]: Set the output format.
          Default is 0.
	  0 - Print output in 3 rows. 
		Row 1: Number_of_rows[\t]Number_of_Columns[\t]Score   
		Row 2: Gene_name_1[b]Gene_name_2[b] ...    
		Row 3: Sample_name_1[b]Sample_name_2[b] ... 

	  1 - Print output in 2 rows (Bibench supported format). 
		Row 1: Row_index_1[b]Row_index_2[b] ...    
		Row 2: Column_index_1[b]Column_index_2[b] ...     


Example: 

1. $./ccsbc -t 0.8 -i ./Synthetic_data_results/Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output_standard.txt -m 50 -p 0

2. $./ccsbc -t 0.8 -i ./Synthetic_data_results/Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output_bibench.txt -m 50 -p 1

See your results at "TempResults"


C. Compilation of CUDA C code
====================================
* Note that a CUDA supported GPU card and CUDA C compile is required CUDA C code  

1. Change your current directory to "CUDA_C".
	$cd CUDA_C

2. Type following in the linux command line 
   	$nvcc ./src/ccs.cu -lm -o ccs_cuda

D. Execute CUDA C code
====================================

Type following in the linux commandprompt

	./ccs_cuda -t [correlation threshold] -i [input file] -o [output file] -m [number_of_base_gene]

Example
	$./ccs_cuda -t 0.8 -i ../Synthetic_data_results/Data/Data_Constant_100_1_bicluster.txt -o ./Output.txt

	$./ccs_cuda -t 0.8 -i ../Synthetic_data_results/Data/Data_Constant_100_1_bicluster.txt -o ./Output.txt -m 50




====================================
Authors
====================================

* Anindya Bhattacharya, <anindyamail123@gmail.com>
* Yan Cui, <ycui2@uthsc.edu>

====================================
Contact
====================================

If you have comments or questions, or if you would like to contribute
to the improvement of CCS, please send us an email at anindyamail123@gmail.com and ycui2@uthsc.edu

====================================
License
====================================

This projected is licensed under the terms of the GNU General Public License v3.0.

