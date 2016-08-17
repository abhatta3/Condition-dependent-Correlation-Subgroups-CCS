A. Compilation of C code
====================================

1. Change your current directory to "src".
	$cd src

2. Type "make" to create executable "ccsbc" in the home directory
	$make


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

Example: $./ccsbc -t 0.8 -i ./Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output.txt -m 50


-p [0/1]: Set the output format.
          Default is 0.
	  0 - Print output in 3 rows. 
		Row 1: Number_of_rows[\t]Number_of_Columns[\t]Score   
		Row 2: Gene_name_1[b]Gene_name_2[b] ...    
		Row 3: Sample_name_1[b]Sample_name_2[b] ... 

	  1 - Print output in 2 rows (Bibench supported format). 
		Row 1: Row_index_1[b]Row_index_2[b] ...    
		Row 2: Column_index_1[b]Column_index_2[b] ...     


Example: $./ccsbc -t 0.8 -i ./Data/Data_Constant_100_1_bicluster.txt -o ./Results/Output.txt -m 50 -p 1




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
	$./ccs_cuda -t 0.8 -i ../Data/Data_Constant_100_1_bicluster.txt -o ./Output.txt

	$./ccs_cuda -t 0.8 -i ../Data/Data_Constant_100_1_bicluster.txt -o ./Output.txt -m 50



