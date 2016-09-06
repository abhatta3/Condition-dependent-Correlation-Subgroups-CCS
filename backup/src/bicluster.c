#include "../include/ccs.h"

void computebicluster(struct gn *gene, int n,int D,double thr,int k,struct bicl *maxbc)
{

	struct bicl tmpbc;
	double tempj,mean_k,mean_i;
	int i,j,l,vl,lc,tc,wid[3];
	char *vect[3];
	int *sxc,c;
	int dif,ta,ii,overlapflag=0;
	double sx, sxx,	r, sy, sxy, syy, *sxxb, *sxb;
	double mintr, tr, tsy, tsyy, tsxy;


  	for (i = 0; i < 3; i++)
  		vect[i] = (char *)calloc(D,sizeof(char));
  
  
	tmpbc.sample = (char *)calloc(D,sizeof(char));
	tmpbc.data = (char *)calloc(n,sizeof(char));

	sxxb = (double *)calloc(n,sizeof(double)); 
	sxb = (double *)calloc(n,sizeof(double));
	sxc = (int *)calloc(n,sizeof(int));
        

	maxbc[k].score=1.0;

       //calculate mean expression for gene k
        mean_k=gene[k].x[D];

  	for (i = k+1; i < n; i++) //pair k,i
    	{ 	
               //calculate mean expression for gene i
		mean_i=gene[i].x[D];

      		wid[0]=0; wid[1]=0; wid[2]=0;      

		for (j = 0; j < D; j++)  
       		{

    	 		if ((gene[k].x[j] - mean_k)>=0 && (gene[i].x[j] - mean_i)>=0) //i and k upregulated : positive correlation
    	 		{
          
          			vect[0][j] = '1';
          			vect[1][j] = '0';
          			vect[2][j] = '0';

                                wid[0]+=1;
         		}
         		else if ((gene[k].x[j] - mean_k)<0 && (gene[i].x[j] - mean_i)<0)  // i and k down regulated : positive correlation
    	 		{
          
          			vect[0][j] = '0';
          			vect[1][j] = '1';
          			vect[2][j] = '0';

                                wid[1]+=1;
         		}
                        else if ((gene[k].x[j] - mean_k)*(gene[i].x[j] - mean_i)<0) //betwenn i and k one is up regulated and the other one is down regulated : negative correlation
                        {

          			vect[0][j] = '0';
          			vect[1][j] = '0';
          			vect[2][j] = '1';
                                wid[2]+=1;

                        } 
  
	       }

      	       for (vl = 0; vl < 3; vl++)
      	       { 
                   	  
                   if(wid[vl]>=min) 
                   {
      			sx = 0;
      			sxx = 0;

      			r = 0;
      			sy = 0;
      			sxy = 0;
      			syy = 0;

      			for (j = 0; j < D; j++) {
				if(vect[vl][j]=='1')
	      	 			sx +=  gene[k].x[j];
		  	}
      		        sx /= wid[vl];
      			for (j = 0; j < D; j++) {
				if(vect[vl][j]=='1')
	      	 			sxx += (sx-gene[k].x[j]) * (sx-gene[k].x[j]);
			}
      			sxx = (double)sqrt(sxx);
           
      			for (j = 0; j < D; j++) {
				if(vect[vl][j]=='1')
	      	 			sy +=  gene[i].x[j];
			}
      			sy /= wid[vl]; 
      
			for (j = 0; j < D; j++)
      			{

				if(vect[vl][j]=='1') {

	        			sxy += (sx - gene[k].x[j]) * (sy - gene[i].x[j]);
        				syy += (sy - gene[i].x[j]) * (sy - gene[i].x[j]);
				}
      			}
      
			syy = (double)sqrt(syy);
   
      			r = sxy/(sxx * syy);

        		if(r<0)
         			r=0-r;

		    }
                    else
                    {
                        r=0.0;
                        continue;
                    }
      			
	
     		    if (r > thr) 
      		    {
          		for (j = 0;j < D; j++)
             			tmpbc.sample[j] = vect[vl][j];
			for (j = 0;j < n; j++)
             			tmpbc.data[j] = '0';
          		tmpbc.data[k] = '1';
          		tmpbc.data[i] = '1';
          		tmpbc.datacount = 2;
          		tmpbc.samplecount = wid[vl];

          		sxc[0]=k;
	  		sxc[1]=i;

          		sxb[0]=sx;
          		sxb[1]=sy;

          		sxxb[0]=sxx;
          		sxxb[1]=syy;
			
          		for (l = 0; l < n; l++)  //bicluster augmentation
          		{

            			if (l != i && l != k)
            			{
              				mintr=1.0; tr = 0; tsy = 0; tsyy = 0;
   			              	for (j = 0; j < D; j++) {
						if(vect[vl][j]=='1')	
	            	  				tsy += gene[l].x[j];
					}
              				tsy /= wid[vl]; 

	       				for (j = 0; j < D; j++){
						if(vect[vl][j]=='1')
	            	   				tsyy += (tsy-gene[l].x[j]) * (tsy-gene[l].x[j]);
              				}
              				tsyy = (double)sqrt(tsyy);

              				for (lc = 0; lc < tmpbc.datacount; lc++)  //correlation with the genes in the bicluster[count]
              				{  
                				tc=sxc[lc]; 
                                                tsxy = 0; 
            					for (j = 0; j < D; j++){
							if(vect[vl][j]=='1')	
	              						tsxy += (sxb[lc]-gene[tc].x[j]) * (tsy-gene[l].x[j]);
            					}
            					tr = tsxy / (sxxb[lc] * tsyy); 
            					if(tr < 0)
               						tr=0-tr;
                				if(tr<=mintr)
                    					mintr=tr;
                 
             				}
                                       
             				if (mintr > thr)
             				{
                 				tmpbc.data[l] = '1';
                  				tmpbc.datacount++;
                 				sxb[lc]=tsy;  
                 				sxxb[lc]=tsyy;
                 				sxc[lc]=l;
             				}         	
          			}
        		}  // end of augmentation

			///////////////////////////Compute number of correlation for other samples////////////////////////

			tc=0;
			dif=0;
			for(j=0;j<n;j++) {

				if(tmpbc.data[j]=='1')
				{
					
	      				sx = 0;
	      				sxx = 0;
		      			for (ii = 0; ii < D; ii++) {
						if(vect[vl][ii]=='0')
			      	 			sx +=  gene[j].x[ii];
                                                        tc++; 
		  			}
                                        if (tc==0)
                                            sx=0;
                                        else
		      		            sx /= tc;
		      			for (ii = 0; ii < D; ii++) {
						if(vect[vl][ii]=='0')
	      			 			sxx += (sx-gene[j].x[ii]) * (sx-gene[j].x[ii]);
					}
		      			sxx = (double)sqrt(sxx);
	
					for(l=j+1;l<n;l++) {
						if(tmpbc.data[l]=='1')
						{
							
				      			r = 0;
				      			sy = 0;
				      			sxy = 0;
				      			syy = 0;
	
      							for (ii = 0; ii < D; ii++) {
								if(vect[vl][ii]=='0')
			      			 			sy +=  gene[l].x[ii];
							}
                                                        if (tc==0)
                                                            sy=0; 
                                                        else
				      			    sy /= tc; 
      
							for (ii = 0; ii < D; ii++)
      							{
								if(vect[vl][ii]=='0') {
				        				sxy += (sx - gene[j].x[ii]) * (sy - gene[l].x[ii]);
			        					syy += (sy - gene[l].x[ii]) * (sy - gene[l].x[ii]);
								}
			      				}

							syy = (double)sqrt(syy);
			      				r = sxy/(sxx * syy);
			        			if(r<0)
			         				r=0-r;
							if(r>thr)
								dif++;
						}
					}
				}
			}
                   
			////////////////////Compute Jaccard score//////////////////////////

			tempj=dif/(tmpbc.datacount*(tmpbc.datacount-1)/2.0);   
                        
                        //////////////////////////Select bicluster candidate as the largest (maxbc[k].datacount<tmpbc.datacount)/////////////////////// 
                        /////////////////////////of all condition dependent (jaccard score <0.01) bicluster for k///////////////////////////////////////   		

                        if(tempj<0.01 && maxbc[k].datacount<tmpbc.datacount && tmpbc.datacount>mingene)
	                {
				maxbc[k].score=tempj;
	 	   		for (j = 0; j < n; j++)  
	 	               		maxbc[k].data[j]=tmpbc.data[j];
 	   			for (j = 0; j < D; j++)  
 	                   		maxbc[k].sample[j]=tmpbc.sample[j];
                        	maxbc[k].datacount=tmpbc.datacount;
                        	maxbc[k].samplecount=tmpbc.samplecount;
				
	                 }

      		} //end of r>thr condition
            } //end of loop for vl

    	}// end of i loop


 	for (i = 0; i < 3; i++)
  		free(vect[i]);
  
 
	free(tmpbc.sample);
	free(tmpbc.data);

	free(sxxb); 
	free(sxb);
	free(sxc);
        


}
