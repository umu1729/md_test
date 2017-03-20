//######################################################
//# this is the code of 'BlockLevel Prefix Sum(Scan)'. #
//######################################################


#pragma once
#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>
#include "common/common_vc.h"


//Sum from smem[0] to smem[blockDim.x-1]. 
//need to syncthreads() before use this device function.(if smem is on SharedMemory.)
__device__ inline float blockLevelSum(float *smem){//only for 2^n (2^n>32)
    //int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int sdx = threadIdx.x;
    int wdx = sdx%32;
    int bDim = blockDim.x;

    //warp level
    float tmp = smem[sdx];
    
    tmp += __shfl_xor(tmp,1);
    tmp += __shfl_xor(tmp,2);
    tmp += __shfl_xor(tmp,4);
    tmp += __shfl_xor(tmp,8);
    tmp += __shfl_xor(tmp,16);
    
    if (wdx == 0) smem[sdx] = tmp;
    
    //blocklevel
    for (int i=5;(1<<(i+1))<=bDim;i++){
        int pos = sdx<<(i+1);
        int wid = 1<<i;
        float tmp = 0.f;
        if ( pos < bDim ){
            tmp = smem[pos]+smem[pos+wid];
        }
        __syncthreads();
        if ( pos < bDim){
            smem[pos] = tmp;
        }
        __syncthreads();
    }
    
    return smem[0];
}

__global__ void blockSum(float* s,float *t){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int sdx = threadIdx.x;
    
    extern __shared__ float smem[];
   
    smem[sdx] = s[idx]; 
    __syncthreads();
    
    float sum = blockLevelSum(smem);
    t[idx] = sum;
    


}

void initRand(){
    srand((unsigned) time(nullptr));
}

void initHBuf(float *buf,int nElem){
    for (int i=0;i<nElem;i++){
        buf[i]=1;
    }
}

int main(int argc,char** argv){
    initRand();
    
    int nElem = 1024*1024;
    size_t nByte = nElem * sizeof(float);
    
    float *h_s,*h_t;
    float *d_s,*d_t;
    
    h_s = (float*)malloc(nByte);
    h_t = (float*)malloc(nByte);
    
    cudaMalloc(&d_s,nByte);
    cudaMalloc(&d_t,nByte);
    
    initHBuf(h_s,nElem);
    cudaMemcpy(d_s,h_s,nByte,cudaMemcpyHostToDevice);
    blockSum<<<nElem/1024,1024,1024*sizeof(float)>>>(d_s,d_t);
    cudaMemcpy(h_t,d_t,nByte,cudaMemcpyDeviceToHost);
    
        CHECK(cudaGetLastError());
    
    for (int i=0;i<10;i++){
        int j = i * 1024;
        printf("%f,%f\n",h_s[j],h_t[j]);
    }
    
    free(h_s);
    free(h_t);
    cudaFree(d_s);
    cudaFree(d_t);
}