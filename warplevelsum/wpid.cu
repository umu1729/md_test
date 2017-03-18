#pragma once
#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>
#include <math.h>
#include "common/common_vc.h"


__global__ void warpLevelSum(int* s,int *t){
    int id = threadIdx.x + blockIdx.x * blockDim.x;
    int sd = threadIdx.x;
    int wd = sd%32;
    
    int tmp,swp;
    tmp = s[id];
    tmp = tmp != __shfl_down(tmp,1) ? 1 : 0;
    
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

void initHBuf(int *buf,int nElem){
    int h = 0;
    for (int i=0;i<nElem;i++){
        int j = 323322 + i;
        h=((j*13)%7+(h*17)%5)%2 ;
        buf[i]=h;
    }
}

int main(int argc,char** argv){
    initRand();
    
    int nElem = 1024*1024;
    size_t nByte = nElem * sizeof(float);
    
    int *h_s,*h_t;
    int *d_s,*d_t;
    
    h_s = (int*)malloc(nByte);
    h_t = (int*)malloc(nByte);
    
    cudaMalloc(&d_s,nByte);
    cudaMalloc(&d_t,nByte);
    
    initHBuf(h_s,nElem);
    cudaMemcpy(d_s,h_s,nByte,cudaMemcpyHostToDevice);
    warpLevelSum<<<nElem/1024,1024>>>(d_s,d_t);
    cudaMemcpy(h_t,d_t,nByte,cudaMemcpyDeviceToHost);
    
    for (int i=0;i<1024;i++){
        printf("%d,%d\n",h_s[i],h_t[i]);
    }
    
    free(h_s);
    free(h_t);
    cudaFree(d_s);
    cudaFree(d_t);
}