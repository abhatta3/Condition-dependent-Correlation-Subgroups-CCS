/**************************************************************************
 * Complete Correlation Subgroups (CCS) 
 * Description: Biclustering has been emerged as a powerful tool for 
                identification of a group of co-expressed genes under a subset 
                of experimental conditions (measurements) present in a gene 
                expression dataset.  In this program we implemented CCS biclustering. 
 * Developer: Dr. Anindya Bhattacharya and Dr. Yan Cui, UTHSC, Memphis, TN, USA
 * Email: anindyamail123@gmail.com; ycui2@uthsc.edu 

Note: The minimum number of genes and the samples per bicluster is 10. 
User can alter the minimum size by changing the values for 'mingene' 
and 'min' defined in "ccs.h" file for minimum number of genes and samples
respectively. 

****************************************************************************
*/

#include "ccs_cu.h"
#include "matrixsize.c"
#include "readgene_cu.c"
#include "bicluster.cu"
#include "print_bicluster_cu.c"


int main(int argc, char *argv[])
{
	FILE   *in,*out;
	struct cgn *gene,*device_gene;
	struct cbicl *bicluster,*device_bicluster;
	char **Hd;
	char *infile,*outfile;	
	int c, errflag,maxbcn=MAXB;
	int i,n,D;
	extern char *optarg;
	float thr;
          
	 clock_t start = clock() ;
	 in= out = NULL;

  	 n = D = 0;
	 thr = errflag = 0;
  
	 while ((c = getopt(argc, argv, "ht:m:i:o:?")) != -1)
  	 {
    		switch(c)
    		{
    			case 'h': // help
      				printUsage();
      				exit(0);
    			case 't': // threshold value
      				thr = atof(optarg);
      				break;
    			case 'm': // maximum number of bicluster search
      				maxbcn = atoi(optarg);
      				break;
    			case 'i': // the input expression file
      				infile = optarg;
      				break;
    			case 'o': // the output file
      				outfile = optarg;
      				break;
    			case ':':       /* -f or -o without operand */
      				printf("Option -%c requires an operand\n", optopt);
      				errflag++;
     				break;
    			case '?':
      				fprintf(stderr,"Unrecognized option: -%c\n", optopt);
      				errflag++;
    		}
 	  }



  	if (thr == 0)
  	{
    		fprintf(stderr,"***** WARNING: Threshold Theta (corr coeff) "
                   "value assumed to be ZERO (0)\n");
  	}


  	if (outfile[0] == '\0')
  	{
    		fprintf(stderr,"***** WARNING: Output file assumed to be STDOUT\n");
    		out = stdout;
  	}
  	else if ((out = fopen(outfile,"w")) == NULL) //write open bicluster file
  	{
    		fprintf(stderr,"***** ERROR: Unable to open Output file %s\n",outfile);
    		errflag++;
  	}

	if ((thr < 0) || (thr > 1))
  	{
    		fprintf(stderr,"***** ERROR: Threshold Theta (corr coeff) "
                   "must be between 0.0-1.0\n");
  	}

  	if (infile[0] == '\0')
  	{
    		fprintf(stderr,"***** ERROR: Input file not defined\n");
    		if (out) fclose(out);
    			errflag++;
  	}
 	else if ((in = fopen(infile,"r")) == NULL)  //open gene file
  	{
    		fprintf(stderr,"***** ERROR: Unable to open Input %s\n", infile);
    		if (out) fclose(out);
    			errflag++;
  	}

  	if (errflag)
  	{
    		printUsage();
    		exit(1);
  	}

	getmatrixsize(in,&n,&D);
	printf("Number of rows=%d\tNumber of columns=%d\n",n,D);


	if(maxbcn>n)
		maxbcn=n;


 	gene = (struct cgn *)calloc(n,sizeof(struct cgn));
 	Hd = (char **)calloc(D+1,sizeof(char *));

  	readgene(infile,gene,Hd,n,D);	

	bicluster = (struct cbicl *)calloc(maxbcn,sizeof(struct cbicl));


	cudaMalloc((void**)&device_gene, n*sizeof(struct cgn));
	cudaMalloc((void**)&device_bicluster, maxbcn*sizeof(struct cbicl));
	cudaMemcpy(device_gene, gene, n*sizeof(struct cgn), cudaMemcpyHostToDevice);


	computebicluster_cu<<<maxbcn,1>>>(device_gene,n,maxbcn,D,thr,device_bicluster);

        cudaMemcpy(bicluster, device_bicluster, maxbcn*sizeof(struct cbicl), cudaMemcpyDeviceToHost); 
   
        printbicluster(out,gene,Hd,n,D,maxbcn,thr,bicluster);


              
 	clock_t end = clock() ;
 	double elapsed_time = (end-start)/(double)CLOCKS_PER_SEC ;
 	printf("Ellapsed time=%lf\n",elapsed_time);

 	fprintf(out,"\n\nEllapsed time=%lf\n",elapsed_time);


  	if (out) fclose(out);


	for (i = 0; i < D+1; i++)   {   
		free(Hd[i]);
	}
	free(Hd);
	free(gene);
	free(bicluster);
	cudaFree(device_gene); 
	cudaFree(device_bicluster);

  	return SUCCESS;
}

void printUsage()
{
printf("\n\t\tUsage: affyAddGene\n"
         "\t\t         -h [display this help message]\n"
         "\t\t         -t threshold theta in a range 0.0 - 1.0\n"
         "\t\t         -m maximum expected biclusters in a range 1 - number_of_rows_in_input_data_matrix\n"
         "\t\t         -i <input microarray expression file (processed data)>\n"
         "\t\t         -o <Output file>\n"
        );
}

