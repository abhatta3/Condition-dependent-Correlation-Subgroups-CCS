Publication link
————————————————
https://www.nature.com/articles/s41598-017-04070-4


*This file is part of the CCS package for biclustering analysis

Name: Condition-dependent Correlation Subgroup (CCS)
-----------------------------------------------------

Introduction to CCS
------------------------
Condition-dependent Correlation Subgroup (CCS) is a biclustering algorithm for comprehensive discovery of functionally coherent biclusters from large-scale gene expression data. For the details of the CCS algorithm, see our paper entitled “A GPU-accelerated algorithm for biclustering analysis and detection of condition-dependent coexpression network modules”. The algorithm is implemented in C. See the steps in the section A and B for compilation and execution of the C code.
The structure of the CCS algorithm is particularly suitable for parallel computing. A CUDA C based GPGPU computing code is included here as a parallel version of the algorithm. You need a programmable GPU card and CUDA C complier for compilation and execution of our code. Follow the steps in C and D.
The performance of CCS was tested on synthetic and real gene expression datasets. We also showed that there is an equivalence between CCS biclusters and condition-dependent co-expression network modules. The related Python codes (Python 2.7 or higher) are available in the "python_utility" directory.
Synthetic and real gene expression datasets, and CCS biclustering results are available in the "Results" directory.


Installation and Execution
--------------------------

A. Compilation of C code
------------------------
1.	Change your current directory to "src". $cd src
2.	Type "make" to create executable "ccsbc" in the home directory $make
3.	Back to parent directory cd ../

B. Execute C code
-----------------
Type the following commands in Linux:
./ccs -t [correlation threshold] -i [input file] -o [output file]
Parameters:
-t [0-1.0]: Specify correlation threshold "theta" between 0 to 1. Recommended value is 0.8. -i [input_data_file] -o [output_data_file]
Example: 
./ccs -t 0.8 -i ./Results/Synthetic_data_results/Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output.txt
Additional parameters:
-m [1 - number of gene/rows in the data matrix]: Set the number of base gene that are to be considered for forming biclusters. Default value is 1000 or maximum number of genes when that is less than 1000.
Example: 
./ccs -t 0.8 -i ./Results/Synthetic_data_results/Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output.txt -m 90
-g [0.0 - 100.0]: Minimum gene set overlap required for merging the overlapped biclusters. Default value is 100.0 for 100% overlap.
Example: 
./ccs -t 0.8 -i ./Results/Synthetic_data_results/Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output.txt -g 50.0
-p [0/1]: Set the output format. Default is 0.
0 - Print output in 3 rows. 

    Row 1: Number_of_rows[\t]Number_of_Columns[\t]Score   
    Row 2: Gene_name_1[b]Gene_name_2[b] ...    
    Row 3: Sample_name_1[b]Sample_name_2[b] ... 

1 - Print output in 2 rows (Bibench supported format). 

    Row 1: Row_index_1[b]Row_index_2[b] ...    
    Row 2: Column_index_1[b]Column_index_2[b] ...     

Example: 
./ccs -t 0.9 -i ./Results/Synthetic_data_results/Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output_standard.txt -m 50 -p 1 -g 100.0

C. Compilation of CUDA C code
------------------------------
*Note that a CUDA supported GPU card and CUDA C compiler is required.
1.	Change your current directory to "CUDA_C". $cd CUDA_C
2.	Type following in the linux command line $nvcc ./src/ccs.cu -lm -o ccs_cuda

D. Execute CUDA C code
-----------------------
Type following commands in the Linux: 
./ccs_cuda -t [correlation threshold] -i [input file] -o [output file]
Example: 
./ccs_cuda -t 0.9 -i ../Synthetic_data_results/Data/Data_Constant_100_1_bicluster.txt -o ./Output.txt -m 50 -p 1 -g 100.0


Authors
-------------------------------------------------
•	Anindya Bhattacharya, anindyamail123@gmail.com
•	Yan Cui, ycui2@uthsc.edu

Contact
-------------
If you have comments or questions, or if you would like to contribute to the further development of CCS, please send us an email at anindyamail123@gmail.com and ycui2@uthsc.edu

License
------------
This projected is licensed under the terms of the GNU General Public License v3.0.
