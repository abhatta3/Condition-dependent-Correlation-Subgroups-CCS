#include "../include/ccs.h"


void printbicluster(FILE *out,struct gn *gene,char **Hd, int n,int D,int maxbcn,double thr,struct bicl *maxbc, int print_type)
{
	int i,j,k,kk,max, maxsize;
        double score;

	for (k=0;k<maxbcn;k++)
  	{
                if (maxbc[k].score>=0.01  || maxbc[k].datacount<mingene || maxbc[k].samplecount<min)
                   continue;

     		for (kk=k+1;kk<maxbcn;kk++) {
	                if (maxbc[kk].score>=0.01 || k==kk || maxbc[kk].datacount<mingene || maxbc[kk].samplecount<min)
	                   continue;
                        score=between_bicluster_correlation(gene,maxbc, k,kk,n,D,thr); 
                        if (score<0.01)
                     	    mergebcl(maxbc, k,kk,n,D,score);
                }

        }  

 
	for (k=0;k<maxbcn;k++)
  	{
               
     		if(maxbc[k].score<0.01  && maxbc[k].datacount>mingene && maxbc[k].samplecount>min)
    		{

                    if(print_type==0)  { 
	    		fprintf(out,"%d\t%d\t%lf\n",maxbc[k].datacount,maxbc[k].samplecount,maxbc[k].score);
 
	    		for (i = 0; i < n; i++)//print genes
    			{
		      		if (maxbc[k].data[i] == '1')  {
		        		fprintf(out,"%s ",gene[i].id);
					      	  
                                } 
    			}
    			fprintf(out,"\n");
    			for (j = 0; j < D; j++)  //print samples
    			{
      				if(maxbc[k].sample[j] == '1')
      					fprintf(out,"%s ",Hd[j+1]);	
    			}
    			fprintf(out,"\n");
                    }
                    else   {
	    		for (i = 0; i < n; i++)//print gene index
    			{
		      		if (maxbc[k].data[i] == '1')  {
		        		fprintf(out,"%d ",gene[i].indx);
                                } 
    			}
    			fprintf(out,"\n");
    			for (j = 0; j < D; j++)  //print samples index
    			{
      				if(maxbc[k].sample[j] == '1')
      					fprintf(out,"%d ",j);	
    			}
    			fprintf(out,"\n\n");


                   }


     		}
	}
}
