#pragma once
#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>
#include "common/common_vc.h"


///
/// This is Pseudo-Naive-Bitonic-Sort;
///

__global__ void bitonicSort_Global(int *dep,int *tar,int i,int j){
    int Pad = blockDim.x * gridDim.x;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int dec = (idx >> (i+1))&0x1; //ascend or decend flag;
    int jdx = idx ^ (1<<j); //target buf
    int rig = (idx >> j)&0x1; //left or right (buffer pointer) flag;
    int sel = dep[idx];
    int opp = dep[jdx];
    bool ima = opp < sel; //AM I MAX;
    bool ins = ((dec ^ rig) == ima) | sel == opp; // is it NOT need to swap memory(idx<-->idy)
    
    //int mi = min(sel,opp);
    //int ma = max(sel,opp);
    //tar[idx] = ((dec ^ rig) == 1)?ma:mi;  // this code can't be extend to hash sorting;
    tar[idx] = ins?sel:opp;
    tar[idx+Pad] = ins?dep[idx+Pad]:dep[jdx+Pad];
}


// 2^j <= blockDim  & blockDim = 2^k
__global__ void bitonicSort_Local(int *dep,int *tar,int i,int jMax){
    extern __shared__ int sbuf [];
    int Pad = blockDim.x * gridDim.x;
    int sPad = blockDim.x;
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int itx = threadIdx.x;
    
    sbuf[itx] = dep[idx];
    sbuf[itx+sPad] = dep[idx+Pad];
    __syncthreads();
    
    
    for (int j=jMax;j>=0;j--){
        int dec = (idx >> (i+1))&0x1;
        int jtx = itx ^ (1<<j);
        int rig = (itx >> j)&0x1;
        int sel = sbuf[itx];
        int opp = sbuf[jtx];
        bool ima = opp < sel;
        bool ins = (dec ^ rig) == ima  | sel == opp;
        __syncthreads();
        sbuf[itx] = ins?sel:opp;
        sbuf[itx+sPad] = ins?sbuf[itx+sPad]:sbuf[jtx+sPad];
        __syncthreads();
    }

    tar[idx] = sbuf[itx];
    tar[idx+Pad] = sbuf[itx + sPad];
}

void SwapIntPointer(int** a,int **b){
    int *tmp = *b;
    *b = *a;
    *a = tmp;
}

/// d_buf , d_work = buffer [Hash:nElem,Key:nElem] , nElem = 2 << nLog2Elem;
/// using Hash,sorting Key
void sort(int*& h_d_buf,int*& h_d_work,int nLog2Elem){
    
    CHECK(cudaGetLastError());
    
    int nLogE = nLog2Elem;
    int nElem = 1 << nLogE;
    
    for(int i=0;i<nLogE;i++){
        int j = i;
        int nLogL = 8;
        int nSizL = 1<<nLogL;
        
        for(;j>= 0;j--){
                bitonicSort_Global<<<nElem/nSizL,nSizL>>>(h_d_buf,h_d_work,i,j);
                SwapIntPointer(&h_d_buf,&h_d_work);
                cudaDeviceSynchronize();
                //printf(":G:");
            
                CHECK(cudaGetLastError());
        }
        
        /*
        for(;j>= nLogL;j--){
                bitonicSort_Global<<<(nElem+nSizL-1)/nSizL,nSizL>>>(h_d_buf,h_d_work,i,j);
                SwapIntPointer(&h_d_buf,&h_d_work);
                cudaDeviceSynchronize();
                printf("GLOBAL");
        }
        
        CHECK(cudaGetLastError());
        
        bitonicSort_Local<<<(nElem+nSizL-1)/nSizL,nSizL,2*nSizL*sizeof(int)>>>(h_d_buf,h_d_work,i,j);
        SwapIntPointer(&h_d_buf,&h_d_work);
        cudaDeviceSynchronize();
        */
        
        //printf("%d\n",j);
        CHECK(cudaGetLastError());
    }
}