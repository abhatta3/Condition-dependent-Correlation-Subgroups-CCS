/**************************************************************************
 * Condition-dependent Correlation Subgroups (CCS) 
 * Compilation: nvcc ccs.cu -o ccs_cuda -lm
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


#include "ccs.h"
#include "matrixsize.c"
#include "readgene.c"
#include "pair_cor.c"
#include "bicluster_pair_score.c"
#include "merge_bicluster.c"
#include "print_bicluster.c"
#define MAXSAMPLE 200 //number of samples in the input datamatrix. Fixed here to make static shared memory in global and device function. They dynamic alternative is commented out.


__device__  struct pair_r comput_r_cuda(char *sample,int wid,int k,int i,int D,float **gene)
{

	float sx, sxx, sy, sxy, syy;
	float sx_n, sxx_n, sy_n, sxy_n, syy_n;
	int j;
        __shared__ float genekj[MAXSAMPLE];
       __shared__ float geneij[MAXSAMPLE];
  
        struct pair_r rval;

        rval.r=0.0;
        rval.n_r=0.0; 

	sx = 0; sxx = 0; sy = 0; sxy = 0; syy = 0;

	sx_n = 0; sxx_n = 0;  sy_n = 0;	sxy_n = 0; syy_n = 0;

	for (j = 0; j < D; j++) {
                genekj[j]=gene[k][j];
		if(sample[j]=='1')
			sx +=  genekj[j];
		else
		        sx_n +=  genekj[j];
	}
	sx /= wid;

	sx_n/=(D-wid);

	for (j = 0; j < D; j++) {
		if(sample[j]=='1')
			sxx += (sx-genekj[j]) * (sx-genekj[j]);
		else
			sxx_n += (sx_n-genekj[j]) * (sx_n-genekj[j]);

	}

	sxx = (float)sqrt(sxx);
	sxx_n = (float)sqrt(sxx_n);

	for (j = 0; j < D; j++) {
                geneij[j]=gene[i][j];
		if(sample[j]=='1')
			sy +=  geneij[j];
		else
			sy_n +=  geneij[j];
	}

	sy /= wid; 

	sy_n /= (D-wid); 

	for (j = 0; j < D; j++)
	{

		if(sample[j]=='1') {

			sxy += (sx - genekj[j]) * (sy - geneij[j]);
			syy += (sy - geneij[j]) * (sy - geneij[j]);
		}
		else {

			sxy_n += (sx_n - genekj[j]) * (sy_n - geneij[j]);
			syy_n += (sy_n - geneij[j]) * (sy_n - geneij[j]);
		}

	}

	syy = (float)sqrt(syy);

	syy_n = (float)sqrt(syy_n);

	rval.r =  fabsf(sxy/(sxx * syy));

	rval.n_r =  fabsf(sxy_n/(sxx_n * syy_n));
        
        return rval;

}



__global__ void computebicluster_cu(float **gene, int n,int maxbcn,int D, float thr, char **maxbc_sample,char **maxbc_data,float *maxbc_score,int *maxbc_datacount, int *maxbc_samplecount, char **tmpbc_sample,char **tmpbc_data)
{

   int k=blockIdx.x*blockDim.x+threadIdx.x;

   if(k<maxbcn) {

                //extern __shared__ char arr[];

		float jcc,mean_k,mean_i;
		int i,j,l,vl,wid,wid_0,wid_1,wid_2,l_i,t_tot,t_dif;
		int dif,tot;
		//char *vect[3];
		__shared__  char vect[3][MAXSAMPLE];
		struct pair_r rval;
                int tmpbc_datacount,tmpbc_samplecount;

                float genekj,geneij;

                //for(i=0;i<3;i++)
                    //vect[i]=&arr[i*D];

		maxbc_score[k]=1.0;
		maxbc_datacount[k]=0;  

 
		//calculate mean expression for gene k

		mean_k=gene[k][D];

		for (i = k+1; i < n; i++) //pair k,i
		{ 	
			//calculate mean expression for gene i
			mean_i=gene[i][D];

			wid_0=0; wid_1=0; wid_2=0;      

			for (j = 0; j < D; j++)  
			{
                                genekj=gene[k][j];
                                geneij=gene[i][j];

				if ((genekj - mean_k)>=0 && (geneij - mean_i)>=0) //i and k upregulated : positive correlation
				{

					vect[0][j] = '1';
					vect[1][j] = '0';
					vect[2][j] = '0';

					wid_0+=1;
				}
				else if ((genekj - mean_k)<0 && (geneij - mean_i)<0)  // i and k down regulated : positive correlation
				{

					vect[0][j] = '0';
					vect[1][j] = '1';
					vect[2][j] = '0';

					wid_1+=1;
				}
				else if ((genekj - mean_k)*(geneij - mean_i)<0) //betwenn i and k one is up regulated and the other one is down regulated : negative correlation
				{

					vect[0][j] = '0';
					vect[1][j] = '0';
					vect[2][j] = '1';
					wid_2+=1;

				} 

			}

			for (vl = 0; vl < 3; vl++)
			{ 
				dif=0; tot=0;
                                if(vl==0)
                                    wid=wid_0; 
                                else if(vl==1)
                                    wid=wid_1; 
                                if(vl==2)
                                    wid=wid_2; 

				  
				if(wid>min) { //minimum samples required to form a bicluster module. Default minimum set to 10 in ccs.h   

					    rval=comput_r_cuda(vect[vl],wid, k, i, D, gene);
				}
				else {

					continue;
				}

				if (rval.r > thr) 
				{
					tot++;      
					if(rval.n_r>thr)
					    dif++;

					for (j = 0;j < D; j++)
						tmpbc_sample[k][j] = vect[vl][j];

					for (j = 0;j < n; j++)
						tmpbc_data[k][j] = '0';

					tmpbc_data[k][k] = '1';
					tmpbc_data[k][i] = '1';
					tmpbc_datacount = 2;
					tmpbc_samplecount = wid;


					for (l = 0; l < n; l++)  { //bicluster augmentation
						if (l != i && l != k) {
		                                        t_tot=0; t_dif=0;
		                                        for(l_i=0;l_i<n;l_i++) {
		                                                    if(tmpbc_data[k][l_i]=='1')  {
		                                                            rval=comput_r_cuda(vect[vl],wid, l, l_i, D, gene);
		                                                    
								            if(rval.r>thr) 
									               t_tot+=1;
		                                                            else {
		                                                                       t_tot=0;
		                                                                       break;
		                                                            }   
		                                                            if(rval.n_r>thr) 
										       t_dif+=1;
		                                                       }  
		                                         }                                                                    


							 if(t_tot>0)  {
		                                    	            tmpbc_data[k][l] = '1';
								    tmpbc_datacount+=1;
		                                                    tot+=t_tot; dif+=t_dif;
		                                          }
						}


					}  // end of augmentation

					// Compute Jaccard score

					if(tot>0)
					    jcc=(float)dif/tot;   
					else
					   jcc=1.0; 

					/*   Select bicluster candidate as the largest (maxbc[k].datacount<tmpbc.datacount) 
					     of all condition dependent (jaccard score <0.01) bicluster for k. Minimum number of gene 
		                             for a bicluster is set at 10. See the mingene at ccs.h                                */

					if(jcc<0.01 && maxbc_datacount[k]<tmpbc_datacount && tmpbc_datacount>mingene)
					{
						maxbc_score[k]=jcc;
						for (j = 0; j < n; j++)  
					       		maxbc_data[k][j]=tmpbc_data[k][j];
						for (j = 0; j < D; j++)  
					   		maxbc_sample[k][j]=tmpbc_sample[k][j];
						maxbc_datacount[k]=tmpbc_datacount;
						maxbc_samplecount[k]=tmpbc_samplecount;

					 }


				}    //end of r>thr condition
			}    //end of loop for vl	

		}  // end of i loop
         //for(i=0;i<3;i++)
              //free(vect[i]);

	}
}




int main(int argc, char *argv[])
{
	FILE   *in,*out;
	struct gn *gene;
	char **Hd;
	char *infile,*outfile;	
	int c, errflag,maxbcn=MAXB,print_type=0;
	int i,n,D;
	extern char *optarg;
	float thr,**device_gene,**device_gene_temp;
        struct bicl *bicluster;
        float overlap=100.0; 

        char **device_bicluster_sample,**device_bicluster_temp_sample,**device_bicluster_data,**device_bicluster_temp_data;
        float *device_bicluster_score,*device_bicluster_temp_score;

        int *device_bicluster_datacount,*device_bicluster_temp_datacount;
        int *device_bicluster_samplecount,*device_bicluster_temp_samplecount;

        char **device_bicluster_sample_tmpbc,**device_bicluster_temp_sample_tmpbc,**device_bicluster_data_tmpbc,**device_bicluster_temp_data_tmpbc;



	clock_t start = clock() ;
	in= out = NULL;

  	n = D = 0;
	thr = errflag = 0;
  
	while ((c = getopt(argc, argv, "ht:m:i:p:o:g:?")) != -1)
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

    			case 'g': // output file format
      				overlap = atof(optarg);
     				break;

    			case 'p': // output file format
      				print_type = atoi(optarg);
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

 	gene = (struct gn *)calloc(n,sizeof(struct gn));
 	Hd = (char **)calloc(D+1,sizeof(char *));

        cudaMalloc((void **)&device_gene, n*sizeof(float *));
        device_gene_temp = (float **)calloc(n,sizeof(float *));

  	for (i = 0; i < n; i++) {
	    gene[i].x = ( float *)calloc(D+1,sizeof( float));
            cudaMalloc( (void **)&device_gene_temp[i], (D+1)*sizeof(float));

	    if (!gene[i].x) {
		    printf("***** Memory ERROR: Can't allocate memory to read input data; Exiting \n");
		    exit(1);
   	    }
         }

  	readgene(infile,gene,Hd,n,D);	

  	for (i = 0; i < n; i++) {
              cudaMemcpy(device_gene_temp[i],gene[i].x, (D+1)*sizeof(float), cudaMemcpyHostToDevice);
        }  	

        cudaMemcpy(device_gene, device_gene_temp, n*sizeof(float *), cudaMemcpyHostToDevice);


	bicluster = (struct bicl *)calloc(maxbcn,sizeof(struct bicl));

        cudaMalloc((void **)&device_bicluster_sample, maxbcn*sizeof(char *));
        device_bicluster_temp_sample = (char **)calloc(maxbcn,sizeof(char *));

        cudaMalloc((void **)&device_bicluster_data, maxbcn*sizeof(char *));
        device_bicluster_temp_data = (char **)calloc(maxbcn,sizeof(char *));

        cudaMalloc((void **)&device_bicluster_score, maxbcn*sizeof(float));
        device_bicluster_temp_score = (float *)calloc(maxbcn,sizeof(float));


        cudaMalloc((void **)&device_bicluster_datacount, maxbcn*sizeof(int));
        device_bicluster_temp_datacount = (int *)calloc(maxbcn,sizeof(int));

        cudaMalloc((void **)&device_bicluster_samplecount, maxbcn*sizeof(int));
        device_bicluster_temp_samplecount = (int *)calloc(maxbcn,sizeof(int));




        cudaMalloc((void **)&device_bicluster_sample_tmpbc, maxbcn*sizeof(char *));
        device_bicluster_temp_sample_tmpbc = (char **)calloc(maxbcn,sizeof(char *));

        cudaMalloc((void **)&device_bicluster_data_tmpbc, maxbcn*sizeof(char *));
        device_bicluster_temp_data_tmpbc = (char **)calloc(maxbcn,sizeof(char *));



        if (!bicluster) {
		    printf("***** Memory ERROR: Can't allocate memory for biclusters; Exiting \n");
		    exit(1);
	}



  	for (i = 0; i < maxbcn; i++)
  	{
		bicluster[i].sample = (char *)calloc(D,sizeof(char));
		bicluster[i].data = (char *)calloc(n,sizeof(char));

		cudaMalloc( (void **)&device_bicluster_temp_sample[i], D*sizeof(char));
		cudaMalloc( (void **)&device_bicluster_temp_data[i], n*sizeof(char));


		cudaMalloc( (void **)&device_bicluster_temp_sample_tmpbc[i], D*sizeof(char));
		cudaMalloc( (void **)&device_bicluster_temp_data_tmpbc[i], n*sizeof(char));

		if (!bicluster[i].sample || !bicluster[i].data) {
			    printf("***** Memory ERROR: Can't allocate memory for %d thbiclusters; Exiting \n",i+1);
			    exit(1);
		}

  	}
        
         
        cudaMemcpy(device_bicluster_sample, device_bicluster_temp_sample, maxbcn*sizeof(char *), cudaMemcpyHostToDevice);

        cudaMemcpy(device_bicluster_data, device_bicluster_temp_data, maxbcn*sizeof(char *), cudaMemcpyHostToDevice);


        cudaMemcpy(device_bicluster_sample_tmpbc, device_bicluster_temp_sample_tmpbc, maxbcn*sizeof(char *), cudaMemcpyHostToDevice);

        cudaMemcpy(device_bicluster_data_tmpbc, device_bicluster_temp_data_tmpbc, maxbcn*sizeof(char *), cudaMemcpyHostToDevice);

        //int sherememsize=(3*D)*(sizeof(char)); 

        //computebicluster_cu<<<maxbcn,1,sherememsize>>>(device_gene,n,maxbcn,D,thr,device_bicluster_sample,device_bicluster_data,device_bicluster_score,device_bicluster_datacount,device_bicluster_samplecount, device_bicluster_sample_tmpbc,device_bicluster_data_tmpbc);


        computebicluster_cu<<<maxbcn,1>>>(device_gene,n,maxbcn,D,thr,device_bicluster_sample,device_bicluster_data,device_bicluster_score,device_bicluster_datacount,device_bicluster_samplecount, device_bicluster_sample_tmpbc,device_bicluster_data_tmpbc);

        cudaMemcpy(device_bicluster_temp_score, device_bicluster_score, maxbcn*sizeof(float), cudaMemcpyDeviceToHost);
        cudaMemcpy(device_bicluster_temp_datacount, device_bicluster_datacount, maxbcn*sizeof(int), cudaMemcpyDeviceToHost);
        cudaMemcpy(device_bicluster_temp_samplecount, device_bicluster_samplecount, maxbcn*sizeof(int), cudaMemcpyDeviceToHost);


        for(i=0; i<maxbcn; i++)
        {
             cudaMemcpy(bicluster[i].sample, device_bicluster_temp_sample[i], D*sizeof(char), cudaMemcpyDeviceToHost);
             cudaMemcpy(bicluster[i].data, device_bicluster_temp_data[i], n*sizeof(char), cudaMemcpyDeviceToHost);

             bicluster[i].score=device_bicluster_temp_score[i];

             bicluster[i].datacount=device_bicluster_temp_datacount[i];
             bicluster[i].samplecount=device_bicluster_temp_samplecount[i];

         }



        printbicluster(out,gene,Hd,n,D,maxbcn,thr,bicluster,print_type,overlap);




  	for (i = 0; i < n; i++) {
                     
		free(gene[i].x);
		cudaFree(device_gene_temp[i]);
        }

	for (i = 0; i < D+1; i++)   {   
		free(Hd[i]);
	}

	free(Hd);
	free(gene);
        free(device_gene_temp);
        cudaFree(device_gene);

   
	for (i = 0; i < maxbcn; i++)
  	{  
		free(bicluster[i].sample);
		free(bicluster[i].data);
		cudaFree(device_bicluster_temp_sample[i]);
		cudaFree(device_bicluster_temp_data[i]);

		cudaFree(device_bicluster_temp_sample_tmpbc[i]);
		cudaFree(device_bicluster_temp_data_tmpbc[i]);

        }


	cudaFree(device_bicluster_sample);
        cudaFree(device_bicluster_data);
        cudaFree(device_bicluster_score);
        cudaFree(device_bicluster_datacount);
        cudaFree(device_bicluster_samplecount);

	cudaFree(device_bicluster_sample_tmpbc);
        cudaFree(device_bicluster_data_tmpbc);



	free(device_bicluster_temp_data);
	free(device_bicluster_temp_sample);
	free(device_bicluster_temp_score);

	free(device_bicluster_temp_data_tmpbc);
	free(device_bicluster_temp_sample_tmpbc);

        free(bicluster);



 	clock_t end = clock() ;
 	float elapsed_time = (end-start)/( float)CLOCKS_PER_SEC ;
 	printf("Ellapsed time= %f\n",elapsed_time);
        if(print_type==0)     
        	fprintf(out,"\n\nEllapsed time= %f\n",elapsed_time);

  	if (out) fclose(out);

   	return SUCCESS;
}

void printUsage()
{
printf("\n\t\tUsage: ./ccs_cuda\n"
         "\t\t         -h [display this help message]\n"
         "\t\t         -t threshold theta in a range 0.0 - 1.0\n"
         "\t\t         -o output file\n"
         "\t\t         -i input microarray expression file (processed data)\n"
         "\t\t         -m optional parameter maximum expected biclusters in a range 1 - number_of_rows_in_input_data_matrix with a default 1000\n"
         "\t\t         -g optional parameter minimum percentage of gene overlap between merged biclusters in a range 0-100 with a default 100 percent for full overlap\n"
         "\t\t         -p optional output file format : 0 for standard format, 1 for BiBench bicluster format and default is 0\n"
        );
}

