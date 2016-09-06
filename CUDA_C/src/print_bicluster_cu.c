#include "ccs_cu.h"


float between_bicluster_correlation(struct cgn *gene,struct cbicl *bicluster, int nw,int old,int n,int D,float thr)
{
 int i,j,l;
 float r[2],sx[2],sxx[2],sy[2],syy[2],sxy[2],score;
 char *g,*x;
 int tc[2],dif,un,overlap=0;


  for (i=0;i<n;i++) {
      if (bicluster[nw].data[i]=='1' && bicluster[old].data[i]=='1') {
           overlap=1; 
           break;
      }
  }

 if (overlap==0)
      return(1.0);  

 x = (char *)calloc(D,sizeof(char));

 g = (char *)calloc(n,sizeof(char));

 
 //////////////////////////////////////////////////////////////

 for (i=0;i<n;i++) {

  if (bicluster[nw].data[i]=='1' || bicluster[old].data[i]=='1') 
      g[i]='1';
  else
      g[i]='0';
 }


 tc[0]=0; tc[1]=0; 
 for (i=0;i<D;i++) {

  if (bicluster[nw].sample[i]=='1' || bicluster[old].sample[i]=='1') {
      x[i] = '1';
      tc[0]++;
  } 
  else {
      x[i] = '0';
      tc[1]++;
  }
 }

 ///////////////////////////Compute number of correlation for other samples////////////////////////

 
 dif=0;
 un=0;
 for(i=0;i<n;i++) {
 	if(g[i]=='1')
	{
	      sx[0] = 0; sx[1]=0;
	      sxx[0] = 0; sxx[1]=0;


	      for (j = 0; j < D; j++) {
			if(x[j]=='1')
			{ 	
                             sx[0] +=  gene[i].x[j];
                             
		  	}
                        else {
                             sx[1] +=  gene[i].x[j];
                             

                        }
              }
              if (tc[0]==0)
                     sx[0]=0;
              else
		     sx[0] /= tc[0];
              if (tc[1]==0)
                     sx[1]=0;
              else
	             sx[1] /= tc[1];
              

              for (j = 0; j < D; j++) {
		  if(x[j]=='1') {
			sxx[0] += (sx[0]-gene[i].x[j]) * (sx[0]-gene[i].x[j]);
                   }
                   else {
			sxx[1] += (sx[1]-gene[i].x[j]) * (sx[1]-gene[i].x[j]);
                   } 
	      }
	      sxx[0] = (float)sqrt(sxx[0]);
	      sxx[1] = (float)sqrt(sxx[1]);

	      for(l=i+1;l<n;l++) {

    		    if(g[l]=='1')
		    {
						
			r [0] = 0; r [1] = 0;
			sy[0] = 0; sy[1] = 0;
			sxy[0] = 0; sxy[1] = 0;
			syy[0] = 0; syy[1]=0;

  			for (j = 0; j < D; j++) {
				if(x[j]=='1') {
   	 			     sy[0] +=  gene[l].x[j];
                                }
                                else {

   	 			     sy[1] +=  gene[l].x[j];

                                }
   
                         }

                         if (tc[0]==0)
                             sy[0]=0; 
                         else
			     sy[0] /= tc[0]; 
                         if (tc[1]==0)
                             sy[1]=0; 
                         else
   		             sy[1] /= tc[1]; 


                         for (j = 0; j < D; j++)
      		 	 {
			           if(x[j]=='1') {
				        sxy[0] += (sx[0] - gene[i].x[j]) * (sy[0] - gene[l].x[j]);
			        	syy[0] += (sy[0] - gene[l].x[j]) * (sy[0] - gene[l].x[j]);
				   }
                                   else {
				        sxy[1] += (sx[1] - gene[i].x[j]) * (sy[1] - gene[l].x[j]);
			        	syy[1] += (sy[1] - gene[l].x[j]) * (sy[1] - gene[l].x[j]);

                                   }     
			  }

			  syy[0] = (float)sqrt(syy[0]);

			  syy[1] = (float)sqrt(syy[1]);

 			  r[0] = sxy[0]/(sxx[0] * syy[0]);

 			  r[1] = sxy[1]/(sxx[1] * syy[1]);

			  if(r[0]<0)
   				r[0]=0-r[0];
					
			  if(r[1]<0)
   				r[1]=0-r[1];


                          if(r[0]>thr && r[1]>thr)
				dif++;
                          if (r[0]>thr || r[1]>thr)
                                un++;

	          }
	     }
       }
  }
                   


/////////////////////////////////////////////////
if (un>0)
  score=dif/un;   
else 
  score=1.0; 

free(x);
free(g);
                        
 return(score);

}




void mergebcl (struct cbicl *bicluster, int nw, int old, int n, int D,float score)
{
        int i;
   	for (i = 0; i < n; i++) 
        {
      		if (bicluster[old].data[i]=='1' && bicluster[nw].data[i]=='0') 
                {
                    bicluster[nw].data[i]='1';
                    bicluster[nw].datacount+=1;
                }
        }

   	for (i = 0; i < D; i++) 
        {
      		if (bicluster[old].sample[i]=='1' &&  bicluster[nw].sample[i]=='0') 
                {
                    bicluster[nw].sample[i]='1';
                    bicluster[nw].samplecount+=1;
                 }
        } 
        bicluster[old].score=1.0; //Delete the merged bicluster
        bicluster[nw].score=score; //update merged bicluster score 



}



void printbicluster(FILE *out,struct cgn *gene,char **Hd, int n,int D,int maxbcn,float thr,struct cbicl *maxbc)
{
	int i,j,k,kk;
        float score;

	for (k=0;k<maxbcn;k++)
  	{
                if (maxbc[k].score>=0.01 || maxbc[k].datacount<mingene || maxbc[k].samplecount<min)
                   continue;

     		for (kk=k+1;kk<maxbcn;kk++) {
	                if (maxbc[kk].score>=0.01 || k==kk  || maxbc[kk].datacount<mingene || maxbc[kk].samplecount<min)
	                   continue;
                        score=between_bicluster_correlation(gene,maxbc, k,kk,n,D,thr); 

                        if (score<0.01)
                     	    mergebcl(maxbc, k,kk,n,D,score);
                }

        }  

 
	for (k=0;k<maxbcn;k++)
  	{
               
     		if(maxbc[k].score<0.01 && maxbc[k].datacount>mingene && maxbc[k].samplecount>min)
    		{

	    		fprintf(out,"%d\t%d\t%f\n",maxbc[k].datacount,maxbc[k].samplecount,maxbc[k].score);
 
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
	}
}


