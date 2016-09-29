#include "ccs_cu.h"

__global__ void computebicluster_cu(struct cgn *gene, int n,int maxbcn,int D, float thr, struct cbicl *maxbc)
{

   uint k=blockIdx.x*blockDim.x+threadIdx.x;

   if(k<maxbcn) {

	float tempj,mean_k,mean_i;
	int i,j,l,vl,lc,tc,wid[3];
	char vect[3][200];
	int sxc[1500];  //max number of element possible in a bicluster
	int dif,ii;
	float sx, sxx,	r, sy, sxy, syy, sxxb[1500], sxb[1500];
	float mintr, tr, tsy, tsyy, tsxy;
	struct cbicl tmpbc;

	maxbc[k].samplecount=0;
    	maxbc[k].datacount=0;

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
      			sxx = sqrt(sxx);
           
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
      
			syy = sqrt(syy);
   
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

          		tmpbc.data[k]= '1';
          		tmpbc.data[i]= '1';
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
                                if (tmpbc.datacount>=999)
                                     break;
  
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
              				tsyy = sqrt(tsyy);

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
                                                tmpbc.datacount+=1;
                 				sxb[tmpbc.datacount]=tsy;  
                 				sxxb[tmpbc.datacount]=tsyy;
                 				sxc[tmpbc.datacount]=l;
                  				
             				}         	
          			}
        		}  // end of augmentation

			///////////////////////////Compute number of correlation for other samples////////////////////////
			
			dif=0;
                        tc=D-wid[vl];
			for(j=0;j<n;j++) {

				if(tmpbc.data[j]=='1')
				{
	      				sx = 0;
	      				sxx = 0;
		      			for (ii = 0; ii < D; ii++) {
						if(vect[vl][ii]=='0')
			      	 			sx +=  gene[j].x[ii];
                                     
		  			}
                                        if (tc==0)
                                            sx=0;
                                        else
		      		            sx /= tc;
		      			for (ii = 0; ii < D; ii++) {
						if(vect[vl][ii]=='0')
	      			 			sxx += (sx-gene[j].x[ii]) * (sx-gene[j].x[ii]);
					}
		      			sxx = sqrt(sxx);
	
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

							syy = sqrt(syy);
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

        
 
   } //end of k
}
