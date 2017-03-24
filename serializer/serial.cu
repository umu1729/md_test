#pragma once
#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>
#include <string.h>

#include "serial.h"


int main(int argc,char** argv){

    int d1s = 128;
    int* d1 = (int*)malloc(d1s*sizeof(d1s));
    
    for (int i=0;i<d1s;i++){
        d1[i] = i;
    }
    
    int d2s = 64;
    float* d2 = (float*)malloc(d2s*sizeof(float));
    
    int d3s = 16;
    double* d3 = (double*) malloc(d3s*sizeof(double));
    
    for (int i=0;i<d3s;i++){
        d3[i] = ((double)i)/12;
    }
    
    
    ChunkIndex chk1 = constructChunk("int1",d1s*sizeof(int),d1);
    ChunkIndex chk2 = constructChunk("flo1",d2s*sizeof(float),d2);
    ChunkIndex* chksA[2] = {&chk1,&chk2};
    ChunkIndex chk3 = groupChunk("g1",chksA,2);
    ChunkIndex chk4 = constructChunk("dou1",d3s*sizeof(double),d3);
    ChunkIndex* chksB[2] = {&chk3,&chk4};
    ChunkIndex chk5 = groupChunk("g2",chksB,2);
    
    updateChunkInfo(chk5);
    FILE * fp= fopen("dat.txt","wb");
    printf("ch:%d\n",(int)sizeof(ChunkHeader));
    saveChunk(fp,chk5);
    saveChunk(fp,chk5);
    fclose(fp);
    
    printf("LOAD\n");
    
    fp = fopen("dat.txt","rb");
    fseek(fp,0,SEEK_END);
    size_t size = ftell(fp);
    fseek(fp,0,SEEK_SET);
    printf("fsize:%d\n",(int)size);
    void *dat = malloc(size);
    fread(dat,1,size,fp);
    printf("%s\n",NameBuf(dat));
    
    void *tmp = NextBuf(dat);
    printf("%s\n",NameBuf(tmp));
    tmp = DescendBuf(tmp);
    tmp = DescendBuf(tmp);
    printf("%s\n",NameBuf(tmp));
    int nElem = SizeBuf(tmp)/sizeof(int);
    int* tmp2 = (int*)DescendBuf(tmp);
    for (int i=0;i<nElem;i++){
        printf("%d,",tmp2[i]);
    }
    
    tmp = DescendBuf(dat);
    tmp = NextBuf(tmp);
    printf("%s\n",NameBuf(tmp));
    double* tmp3 = (double*)DescendBuf(tmp);
    nElem = SizeBuf(tmp)/sizeof(double);
    for (int i=0;i<nElem;i++){
        printf("%lf,",tmp3[i]);
    }
    
    
    
    return 0;
}