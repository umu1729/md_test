//####################################################
//# this is the code of 'WarpLevel Prefix Sum(Scan)'.#
//####################################################


#pragma once
#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>
#include "common/common_vc.h"


__global__ void warpLevelSum(float* s,float *t){
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    int sd = threadIdx.x;
    int wd = sd%32;
    
    float tmp,swp;
    tmp = s[id];
    
    
    for (int i=0;i<5;i++){
        int p2 = 1<<i;
        int pp2= p2<<1;
        swp = __shfl_xor(tmp,p2);
        tmp +=( wd%pp2 == pp2 -1 )? swp : 0;
    }
    tmp = wd%32 == 31 ? 0 : tmp;
    
    for (int i=4;i>=0;i--){
        int p2 = 1<<i;
        int pp2= p2<<1;
        swp = __shfl_xor(tmp,p2);
        tmp = wd%pp2 == pp2 - p2 - 1 ? swp : tmp;
        tmp += wd%pp2 == pp2 - 1 ? swp : 0;
    }
    
    t[id] = tmp;
}

void initRand(){
    srand((unsigned) time(nullptr));
}

void initHBuf(float *buf,int nElem){
    for (int i=0;i<nElem;i++){
        buf[i]=i+1;
    }
}

int main(int argc,char** argv){
    initRand();
    
    int nElem = 32;
    size_t nByte = nElem * sizeof(float);
    
    float *h_s,*h_t;
    float *d_s,*d_t;
    
    h_s = (float*)malloc(nByte);
    h_t = (float*)malloc(nByte);
    
    cudaMalloc(&d_s,nByte);
    cudaMalloc(&d_t,nByte);
    
    initHBuf(h_s,nElem);
    cudaMemcpy(d_s,h_s,nByte,cudaMemcpyHostToDevice);
    warpLevelSum<<<1,nElem>>>(d_s,d_t);
    cudaMemcpy(h_t,d_t,nByte,cudaMemcpyDeviceToHost);
    
    for (int i=0;i<nElem;i++){
        printf("%f,%f\n",h_s[i],h_t[i]);
    }
    
    free(h_s);
    free(h_t);
    cudaFree(d_s);
    cudaFree(d_t);
}