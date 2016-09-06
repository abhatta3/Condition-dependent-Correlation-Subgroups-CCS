#include "ccs_cu.h"


void readgene(char * infile, struct cgn *gene,char **Hd, int n,int D)
{
	FILE *in;
	int    i,j,l;
	float tr,average;
	char *line;

	line=(char *)malloc(100000*sizeof(char));

	in = fopen(infile,"r");
  
	for (j = 0; j < D+1; j++)   {   
		fscanf(in,"%s",line);
		Hd[j]=(char *)calloc(strlen(line),sizeof(char));
		strcpy(Hd[j],line);
	}

	
	for (i = 0; i < n; i++)	{ 
		fscanf(in,"%s",gene[i].id);
		average=0;
		for (j = 0; j < D; j++)	{
			fscanf(in,"%f",&gene[i].x[j]);
			average+=gene[i].x[j];   
		}
		gene[i].x[D]=average/D;
		average=0;
		for (j = 0; j < D; j++)
			average += (gene[i].x[D]-gene[i].x[j]) * (gene[i].x[D]-gene[i].x[j]);
			gene[i].x[D+1]=average/D;
	}


	for (i = 0; i < n-1; i++)	{  
		for (j = i+1; j < n; j++) {
			if(gene[i].x[D+1]<gene[j].x[D+1]) {
				strcpy(line,gene[j].id);
				strcpy(gene[j].id,gene[i].id);
				strcpy(gene[i].id,line);

				for (l = 0; l < D+2; l++) {  
					tr=gene[i].x[l];
					gene[i].x[l]=gene[j].x[l];
					gene[j].x[l]=tr;
				}
			}
		}
	}



	free(line);
	if (in) fclose(in);
}

