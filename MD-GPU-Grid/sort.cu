#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>


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

__global__ void bitonicSort_Local1(int *dep,int *tar,int i,int j){
    extern __shared__ int sbuf [];
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int itx = threadIdx.x;
    
    sbuf[itx] = dep[idx];
    __syncthreads();
}

void initRand(){
    srand((unsigned) time(nullptr));
}

void initializeBuffer(int *buf,int nElem){
    for (int i=0;i<nElem;i++){
        buf[i]=(int)( ((double)rand()/(RAND_MAX+1))*1000 );
        printf("%03d,",buf[i]);
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
    
    int nLogE = 8;
    int nElem = 1<<nLogE;
    int nBytes = nElem * sizeof(int);
    
    //CPU BUF
    int *h_buf;
    h_buf = (int*)malloc(nBytes);
    
    //GPU BUF
    int *d_buf,*d_bufs;
    cudaMalloc(&d_buf ,nBytes);
    cudaMalloc(&d_bufs,nBytes);
    
    //DATA INIT
    
    initializeBuffer(h_buf,nElem);
    
    //START
    cudaMemcpy(d_buf,h_buf,nBytes,cudaMemcpyHostToDevice);
    
    for(int i=0;i<nLogE;i++){
        for(int j=i;j>=0;j--){
                bitonicSort_Global<<<nElem/32,32>>>(d_buf,d_bufs,i,j);
                SwapIntPointer(&d_buf,&d_bufs);
                cudaDeviceSynchronize();
        }
    }


    
    cudaMemcpy(h_buf,d_buf,nBytes,cudaMemcpyDeviceToHost);
    for (int i=0;i<nElem;i++){
        printf("%03d,",h_buf[i]);
    }printf("\n");
    
    


    return 0;
}