#ifndef _CCS_H_
#define _CCS_H_

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <unistd.h>
#include <ctype.h>
#include <sys/types.h>
#include <math.h>
#include<time.h>

#define SUCCESS                  0
#define TRUE                     1
#define FALSE                    0
#define FAIL                    -1

#define mingene 10    //minimum number of genes in a bicluster 
#define min 10        //minimum number of samples in a bicluster
#define MAXB 1000     //default number of base gene (outer loop) in the decreasing order of SD to consider for forming the biclusters. Overide if -m number>1000 is used 


#define BUFFER_MAXSIZE 1048576


struct gn{
  char *id;
  int indx;
  double *x;
};

struct bicl{
  char *sample,*data;
  int samplecount,datacount; 
  double score;
};


struct pair_r{
double r,n_r;
};

void printUsage();//print the usage of this code for -h option 
void getmatrixsize(FILE *,int *,int *); //read the input datamatrix and find the number of rows n and number of columns D.
void readgene(char *,struct gn *,char **,int,int); //read the input data matrix and put the expression values in a structure with one additional last column that holds the averages. 
struct pair_r comput_r(char *,int,int,int,int,struct gn *); //compute correlation between pair of genes for a given set of samples and samples not in the given set.
void computebicluster(struct gn *,int,int,double,int,struct bicl *,struct bicl,char **); //compute a condition dependent bicluster module
void printbicluster(FILE *,struct gn *,char **,int,int,int,double,struct bicl *,int); // check the bicluster module and call merge for each pair of bicluster module to see if they can be merged. then print the biclusters that have number of rows >mingene and samples >min  
double between_bicluster_correlation(struct gn *,struct bicl *, int,int,int,int,double); //called from printbicluster to get a condition dependent score for merging a pair of bicluster module.
void mergebcl(struct bicl *, int,int,int,int,double); //merge a pair of bicluster module if the sccore from "between_bicluster_correlation" is still satisfy the level set for condition dependency

#endif  // _CCS_H_
