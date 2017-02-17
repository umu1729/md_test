#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>
#include "common/common_vc.h"


///
/// This is Naive-Bitonic-Sort;
///

__global__ void bitonicSort_Global(int *dep,int *tar,int i,int j){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int dec = (idx >> (i+1))&0x1; //ascend or decend flag;
    int jdx = idx ^ (1<<j); //target buf
    int rig = (idx >> j)&0x1; //left or right (buffer pointer) flag;
    int sel = dep[idx];
    int opp = dep[jdx];
    bool ima = opp < sel; //AM I MAX;
    bool ins = ((dec ^ rig) == ima); // is it NOT need to swap memory(idx<-->idy)
    
    //int mi = min(sel,opp);
    //int ma = max(sel,opp);
    //tar[idx] = ((dec ^ rig) == 1)?ma:mi;  // this code can't be extend to hash sorting;
    tar[idx] = ins?sel:opp;
}


// 2^j <= blockDim  & blockDim = 2^k
__global__ void bitonicSort_Local(int *dep,int *tar,int i,int jMax){
    extern __shared__ int sbuf [];
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int itx = threadIdx.x;
    
    sbuf[itx] = dep[idx];
    __syncthreads();
    
    
    for (int j=jMax;j>=0;j--){
        int dec = (idx >> (i+1))&0x1;
        int jtx = itx ^ (1<<j);
        int rig = (itx >> j)&0x1;
        int sel = sbuf[itx];
        int opp = sbuf[jtx];
        bool ima = opp < sel;
        bool ins = (dec ^ rig) == ima;
        __syncthreads();
        sbuf[itx] = ins?sel:opp;
        __syncthreads();
    }

    tar[idx] = sbuf[itx];
}

void initRand(){
    srand((unsigned) time(nullptr));
}

void initializeBuffer(int *buf,int nElem){
    for (int i=0;i<nElem;i++){
        buf[i]=(int)( ((double)rand()/(RAND_MAX+1))*1000 );
        //printf("%03d,",buf[i]);
    }
    printf("\n");
}

void SwapIntPointer(int** a,int **b){
    int *tmp = *b;
    *b = *a;
    *a = tmp;
}


int main(int argc,char** argv){

    initRand();
    
    int nLogE = 17;
    int nElem = 1<<nLogE;
    size_t nBytes = nElem * sizeof(int);
    printf("nElem=%d\n",nBytes);
    
    //CPU BUF
    int *h_buf;
    h_buf = (int*)malloc(nBytes);
    
    //GPU BUF
    int *d_buf,*d_bufs;
    CHECK(cudaMalloc(&d_buf ,nBytes));
    CHECK(cudaMalloc(&d_bufs,nBytes));
    
    //DATA INIT
    
    initializeBuffer(h_buf,nElem);
    
    //START
    CHECK(cudaMemcpy(d_buf,h_buf,nBytes,cudaMemcpyHostToDevice));
    
    for(int i=0;i<nLogE;i++){
        int j = i;
        int nLogL = 8;
        int nSizL = 1<<nLogL;
        for(;j>= nLogL;j--){
                bitonicSort_Global<<<nElem/256,256>>>(d_buf,d_bufs,i,j);
                SwapIntPointer(&d_buf,&d_bufs);
                cudaDeviceSynchronize();
        }
        
        CHECK(cudaGetLastError());
        
        bitonicSort_Local<<<(nElem+nSizL-1)/nSizL,nSizL,nSizL*sizeof(int)>>>(d_buf,d_bufs,i,j);
        SwapIntPointer(&d_buf,&d_bufs);
        cudaDeviceSynchronize();
        
        //printf("%d\n",j);
        CHECK(cudaGetLastError());
    }

    
    
    cudaMemcpy(h_buf,d_buf,nBytes,cudaMemcpyDeviceToHost);
    //for (int i=0;i<nElem;i++){
    //    printf("%03d,",h_buf[i]);
    //}printf("\n");
    
    


    return 0;
}