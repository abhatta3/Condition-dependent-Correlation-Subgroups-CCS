#include "ccs.h"


void computebicluster(struct gn *gene, int n,int D,double thr,int k,struct bicl *maxbc, struct bicl tmpbc,char **vect)
{
	double jcc,mean_k,mean_i;
	int i,j,l,vl,lc,wid[3],l_i,t_tot,t_dif;
	int dif,tot;

	struct pair_r rval;


	maxbc[k].score=1.0;
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
			dif=0; tot=0;
			  
			if(wid[vl]>min) { //minimum samples required to form a bicluster module. Default minimum set to 10 in ccs.h   

				    rval=comput_r(vect[vl],wid[vl], k, i, D, gene);
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
					tmpbc.sample[j] = vect[vl][j];

				for (j = 0;j < n; j++)
					tmpbc.data[j] = '0';

				tmpbc.data[k] = '1';
				tmpbc.data[i] = '1';
				tmpbc.datacount = 2;
				tmpbc.samplecount = wid[vl];


				for (l = 0; l < n; l++)  { //bicluster augmentation
					if (l != i && l != k) {
                                                t_tot=0; t_dif=0;
                                                for(l_i=0;l_i<n;l_i++) {
                                                            if(tmpbc.data[l_i]=='1')  {
                                                                    rval=comput_r(vect[vl],wid[vl], l, l_i, D, gene);
                                                            
						                    if(rval.r>thr && rval.r>rval.n_r) 
							                       t_tot+=1;
                                                                    else {
                                                                               t_tot=0;
                                                                               break;
                                                                    }   
                                                                    if(rval.n_r>thr) 
								               t_dif+=1;
                                                               }  
                                                 }                                                                    


						 if(t_tot>0 && t_tot>t_dif)  {
                                            	            tmpbc.data[l] = '1';
						            tmpbc.datacount++;
                                                            tot+=t_tot; dif+=t_dif;
                                                  }
					}


				}  // end of augmentation

				// Compute Jaccard score

				if(tot>0)
				    jcc=dif/tot;   
				else
				   jcc=1.0; 

				/*   Select bicluster candidate as the largest (maxbc[k].datacount<tmpbc.datacount) 
				     of all condition dependent (jaccard score <0.01) bicluster for k. Minimum number of gene 
                                     for a bicluster is set at 10. See the mingene at ccs.h                                */

				if(jcc<0.01 && maxbc[k].datacount<tmpbc.datacount && tmpbc.datacount>mingene)
				{
					maxbc[k].score=jcc;
					for (j = 0; j < n; j++)  
				       		maxbc[k].data[j]=tmpbc.data[j];
					for (j = 0; j < D; j++)  
				   		maxbc[k].sample[j]=tmpbc.sample[j];
					maxbc[k].datacount=tmpbc.datacount;
					maxbc[k].samplecount=tmpbc.samplecount;

				 }

			}    //end of r>thr condition
		}    //end of loop for vl	

	}  // end of i loop

}
