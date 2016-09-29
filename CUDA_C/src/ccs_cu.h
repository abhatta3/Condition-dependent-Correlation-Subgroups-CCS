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
#define MAXB 1000      //maximum number of biclusters


//For CUDA 
struct cgn{
  char id[100];
  float x[200]; //Maximum 200 samples. If input data is larger than 200 make change this number accordingly 
};

//For CUDA 
struct cbicl{
  char sample[200];
  char data[20000];       //maximum 10,000 genes
  int samplecount,datacount;
  float score;

};


void printUsage();
void getmatrixsize(FILE *,int *,int *);
void readgene(char *,struct cgn *,char **,int,int);
void printbicluster(FILE *,struct cgn *,char **,int,int,int,float,struct cbicl *,int);
void mergebcl(struct cbicl *,int,int,int,int);


#endif  // _CCS_H_
