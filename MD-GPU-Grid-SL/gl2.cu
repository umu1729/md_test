#pragma once
#include <stdio.h>
#include <time.h>
#include "cuHeader.h"
#include <math.h>
#include <cuda_runtime.h>
#include "sortlib.h"

void initRand(){
    time_t t;
    srand((unsigned) time(&t));
}

/// using malloced nElem*3*sizeof(float) buffer,write uniform particles postiion.
/// and Mold buffer pointer to V3Buf
/// nElem must be 4(k)^3 !!! now Implemention is nElem == 256
V3Buf CreateUniformParticles(float* buf,float rho,int nElem,float* p_length){
    int nloop = (int)powf((float)(nElem/4)+0.1f,1.0f/3);
    printf("%d\n",nloop);
    
    float length = powf(nElem/rho,1.0f/3);
    *p_length = length;
    float a = length/nloop;
    float ah = a/2;
    float px[4] = {0,0,ah,ah};
    float py[4] = {0,ah,ah,0};
    float pz[4] = {0,ah,0,ah};
    
    float *h_x,*h_y,*h_z;
    h_x = buf;
    h_y = h_x + nElem;
    h_z = h_y + nElem;
    for (int i=0; i<nElem;i++){
        h_x[i]=0;
        h_y[i]=0;
        h_z[i]=0;
    }
    int i=0;
    for (int ix = 0;ix<nloop;ix++){
        for (int iy = 0;iy<nloop;iy++){
            for(int iz=0;iz<nloop;iz++){
                for (int ia=0;ia<4;ia++){
                    h_x[i] = px[ia] + ix * a;
                    h_y[i] = py[ia] + iy * a;
                    h_z[i] = pz[ia] + iz * a;
                    i++;
                }
            }
        }
    }
    V3Buf v3Buf = {h_x,h_y,h_z,nElem};
    return v3Buf;
}

/// using malloced nElem*3*sizeof(float) buffer,write random particle position.
/// and Mold buffer pointer to V3Buf
/// nElem must be integer.
V3Buf CreateRandomVelocity(float* buf,int nElem){
    float *h_x,*h_y,*h_z;
    h_x = buf;
    h_y = h_x + nElem;
    h_z = h_y + nElem;


    for(int i=0;i<nElem;i++){
        float x = ((float)rand() / ((float)RAND_MAX + 1)) *2 -1;
        float y = ((float)rand() / ((float)RAND_MAX + 1)) *2 -1;
        float z = ((float)rand() / ((float)RAND_MAX + 1)) *2 -1;
        h_x[i]=x;
        h_y[i]=y;
        h_z[i]=z;
    }
    
    V3Buf v3Buf = {h_x,h_y,h_z,nElem};
    return v3Buf;
}

V3Buf CalculateForce(float *force,float *pos,int nElem,float length,double *potential){
    float *h_fx,*h_fy,*h_fz,*h_px,*h_py,*h_pz;
    h_fx = force;
    h_fy = h_fx + nElem;
    h_fz = h_fy + nElem;
    h_px = pos;
    h_py = h_px + nElem;
    h_pz = h_py + nElem;
    for (int i=0;i<nElem*3;i++){//SYOKIKA
        force[i]=0.0f;
    }

    *potential = 0.0;
    
    float eps = 1.0f;
    float sigma = 1.0f;
    float ce12 = 4.0f*eps*powf(sigma,12);
    float ce06 = 4.0f*eps*powf(sigma, 6);
    float cf12 = ce12 * 12.0f;
    float cf06 = ce06 *  6.0f;
    for (int j=0;j<nElem;j++){
        for (int i=0;i<j;i++){
            float dx,dy,dz,r2,r2i,r06i,r12i,fc,fx,fy,fz;
            dx = h_px[i]-h_px[j];
            dy = h_py[i]-h_py[j];
            dz = h_pz[i]-h_pz[j];
            if(dx<-length/2) dx+=length;
            if(dx> length/2) dx-=length;
            if(dy<-length/2) dy+=length;
            if(dy> length/2) dy-=length;
            if(dz<-length/2) dz+=length;
            if(dz> length/2) dz-=length;
            
            r2 = dx*dx+dy*dy+dz*dz;
            if (r2 > 4*4)continue;
            r2i = 1.0f/r2;
            r06i = r2i * r2i * r2i;
            r12i = r06i * r06i;
            *potential += (ce12 * r12i - ce06 * r06i);
            fc = (cf12*r12i-cf06*r06i)*r2i;
            fx = fc * dx;
            fy = fc * dy;
            fz = fc * dz;
            h_fx[i]+=fx;
            h_fy[i]+=fy;
            h_fz[i]+=fz;
            h_fx[j]-=fx;
            h_fy[j]-=fy;
            h_fz[j]-=fz;
        }
    }
    
    
    V3Buf v3Buf = {h_fx,h_fy,h_fz,nElem};
    return v3Buf;
}

double CalculateHamiltonian(float* pos,float* vel,int nElem,double potential){
    float *h_px,*h_py,*h_pz,*h_vx,*h_vy,*h_vz;
    h_px = pos;
    h_py = h_px + nElem;
    h_pz = h_py + nElem;
    h_vx = vel;
    h_vy = h_vx + nElem;
    h_vz = h_vy + nElem;
    
    double energy = 0.0;
    for(int i=0;i<1;i++){
        float px=h_px[i],py=h_py[i],pz=h_pz[i];
        float vx=h_vx[i],vy=h_vy[i],vz=h_vz[i];
        float r = sqrtf(px*px+py*py+pz*pz);
        float v2= vx*vx+vy*vy+vz*vz;
        energy += (double)(v2/2);
    }
    energy += potential;
    //printf("%lf\n",energy);
    
    return energy;
}

void SwapFloatPointer(float** a,float **b){
    float *tmp = *b;
    *b = *a;
    *a = tmp;
}

__global__ void CalculateForce_GPUNaive(float *force,float *pos,int nElem,float length){

    int idx = threadIdx.x + blockIdx.x * blockDim.x; //1dim grid ,1im block

    float *d_fx,*d_fy,*d_fz,*d_px,*d_py,*d_pz;
    d_fx = force;
    d_fy = d_fx + nElem;
    d_fz = d_fy + nElem;
    d_px = pos;
    d_py = d_px + nElem;
    d_pz = d_py + nElem;
    
    //SYOKIKA
    float t_fx=0.0f,t_fy=0.0f,t_fz=0.0f;
    
    //*potential = 0.0; must implement to calculate Hamiltoniam...
    
    //Is it better to load this constants from ConstantMemory?
    float eps = 1.0f;
    float sigma = 1.0f;
    float ce12 = 4.0f*eps*powf(sigma,12);
    float ce06 = 4.0f*eps*powf(sigma, 6);
    float cf12 = ce12 * 12.0f;
    float cf06 = ce06 *  6.0f;
    
    int i = idx;
    for (int j=0;j<nElem;j++){// j->i
        if (i==j)continue;
        float dx,dy,dz,r2,r2i,r06i,r12i,fc,fx,fy,fz;
        dx = d_px[i]-d_px[j];
        dy = d_py[i]-d_py[j];
        dz = d_pz[i]-d_pz[j];
        if(dx<-length/2) dx+=length;
        if(dx> length/2) dx-=length;
        if(dy<-length/2) dy+=length;
        if(dy> length/2) dy-=length;
        if(dz<-length/2) dz+=length;
        if(dz> length/2) dz-=length;
        
        r2 = dx*dx+dy*dy+dz*dz;
        //if (r2 > 4*4)continue; //Cut force with far from cut radius;this may be MEANINGLESS in GPU.
        r2i = 1.0f/r2;
        r06i = r2i * r2i * r2i;
        r12i = r06i * r06i;
        //*potential += (ce12 * r12i - ce06 * r06i);
        fc = (cf12*r12i-cf06*r06i)*r2i;
        fx = fc * dx;
        fy = fc * dy;
        fz = fc * dz;
        t_fx+=fx;
        t_fy+=fy;
        t_fz+=fz;
    }
    d_fx[i]=t_fx;
    d_fy[i]=t_fy;
    d_fz[i]=t_fz;
}

void CalculateForce_UseGPU_Naive(float* h_p,float *h_f,float* d_p,float *d_f,int nElem,float length){
    int nBytes = nElem * sizeof(float);
    cudaMemcpy(d_p,h_p,nBytes*3,cudaMemcpyHostToDevice);
    dim3 block(256);
    dim3 grid((nElem+block.x-1)/block.x);
    CalculateForce_GPUNaive<<<grid,block>>>(d_f,d_p,nElem,length);
    cudaMemcpy(h_f,d_f,nBytes*3,cudaMemcpyDeviceToHost);
}

#define HSF 6//hash size = 1<<HSF

__host__ __device__ inline int getHashPart(float v,int gC,float hL){
    return (int)floorf(v/hL);
}

__host__ __device__ inline int HashParts2Hash(int Hx,int Hy,int Hz){
    return Hx + Hy*(1<<HSF) + Hz*(1<<(2*HSF)); //this implemention uses gC < (1<<HSF);
}

__host__ __device__ inline int getHash(float x,float y,float z,int gC,float hL){
    // gC is grid Count(1dim) 
    // hL is grid Length(one grid,1dim)
    int Hx = getHashPart(x,gC,hL); //!!! if x == hL * gC ,it cause error!!!
    int Hy = getHashPart(y,gC,hL);
    int Hz = getHashPart(z,gC,hL);
    return HashParts2Hash(Hx,Hy,Hz);
}

__host__ __device__ inline void recoverHash(int hash,int &x,int &y,int &z){
    //Not Implement;
}

__global__ void HashFromPosition(float *d_p,int *d_hash,int nElem,int gC,float hL){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    float x = d_p[idx];
    float y = d_p[nElem+idx];
    float z = d_p[2*nElem+idx];
    int h = getHash(x,y,z,gC,hL);
    d_hash[idx] = h;
    d_hash[idx+nElem] = idx; //it's key;
}

__global__ void GenerateHRRM(int *d_hash,int *d_HRRM,int nElem){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int nHash = 1<<(HSF*3);
    
    int ths = d_hash[idx];
    if (idx==0){
        d_HRRM[ths] = 0;
    }
    if(idx+1<nElem){
        int nxt = d_hash[idx+1];
        if (ths != nxt){
            d_HRRM[ths+nHash] = idx+1;
            d_HRRM[nxt]       = idx+1;
        }
    }else{//idx = nElem-1
        d_HRRM[ths+nHash] = nElem;
    }
}

__global__ void CalculateForce_GPUSort(float *force,float *pos,int* hrrm,int *hash,int nElem,float length,int gC,float hL){
    int idx = threadIdx.x + blockIdx.x * blockDim.x; //1dim grid ,1im block
    float *d_fx,*d_fy,*d_fz,*d_px,*d_py,*d_pz;
    d_fx = force;
    d_fy = d_fx + nElem;
    d_fz = d_fy + nElem;
    d_px = pos;
    d_py = d_px + nElem;
    d_pz = d_py + nElem;
    
    //SYOKIKA
    float t_fx=0.0f,t_fy=0.0f,t_fz=0.0f;
    
    //*potential = 0.0; must implement to calculate Hamiltoniam...
    
    //Is it better to load this constants from ConstantMemory?
    float eps = 1.0f;
    float sigma = 1.0f;
    float ce12 = 4.0f*eps*powf(sigma,12);
    float ce06 = 4.0f*eps*powf(sigma, 6);
    float cf12 = ce12 * 12.0f;
    float cf06 = ce06 *  6.0f;
    
    int i = idx;
    float px = d_px[i];
    float py = d_py[i];
    float pz = d_pz[i];
    int Hx = getHashPart(px,gC,hL);
    int Hy = getHashPart(py,gC,hL);
    int Hz = getHashPart(pz,gC,hL);
    int nHash = 1<<(HSF*3);
    
    for (int c=0;c<27;c++){
        int dHx = c%3-1;
        int dHy = (c/3)%3-1;
        int dHz = (c/3/3)%3-1;
        int pHx = (Hx+gC+dHx)%gC;
        int pHy = (Hy+gC+dHy)%gC;
        int pHz = (Hz+gC+dHz)%gC;
        int h = HashParts2Hash(pHx,pHy,pHz);
        int start = hrrm[h];
        int end = hrrm[h+nHash];
        if(!(start<=end & end-start<1000000)){printf("D");};
        
        for (int k=start;k<end;k++){
            int j = hash[k+nElem];
            if (hash[k]!=h)printf("H");
            if (i==j)continue;
            float dx,dy,dz,r2,r2i,r06i,r12i,fc,fx,fy,fz;
            float qx,qy,qz;
            qx = d_px[j];
            qy = d_py[j];
            qz = d_pz[j];
            dx = px-qx;
            dy = py-qy;
            dz = pz-qz;
            if(dx<-length/2) dx+=length;
            if(dx> length/2) dx-=length;
            if(dy<-length/2) dy+=length;
            if(dy> length/2) dy-=length;
            if(dz<-length/2) dz+=length;
            if(dz> length/2) dz-=length;
            if (c==13 & dx>hL*3)printf("G");
            if (getHash(qx,qy,qz,gC,hL)!=h)printf("E");
            r2 = dx*dx+dy*dy+dz*dz;
            r2i = 1.0f/r2;
            r06i = r2i * r2i * r2i;
            r12i = r06i * r06i;
            //*potential += (ce12 * r12i - ce06 * r06i);
            fc = (cf12*r12i-cf06*r06i)*r2i;
            fx = fc * dx;
            fy = fc * dy;
            fz = fc * dz;
            t_fx+=fx;
            t_fy+=fy;
            t_fz+=fz;
        }
    }
    
    d_fx[i]=t_fx;
    d_fy[i]=t_fy;
    d_fz[i]=t_fz;
}

typedef struct{
    int *d_HashKey; //nElem * 2 (pair of hash and key);
    int *d_HashKeyWork;
    int *d_HRRM;//?
    int *h_HashKeyTest;//TestBuf
    int *h_HRRMTest;
    float *h_nElemF;
} WorkBufList;


void CalculateForce_UseGPU(float* h_p,float *h_f,float* d_p,float *d_f,int nElem,float length,WorkBufList wbl){
    float cutRadius = 4.0f;
    int gC = floor(length/cutRadius); //grid Size(1dim)
    float hL = length/gC; //cutRadius
    int nHash = 1<<(HSF*3);
    int nBytes = nElem * sizeof(float);
    
    cudaMemcpy(d_p,h_p,nBytes*3,cudaMemcpyHostToDevice);
    
    
    dim3 block(256);
    dim3 grid(nElem/block.x);
    //printf("G<%d>",grid.x);
    
    CHECK(cudaGetLastError());

    //Kernel:Hashing (Generate Hash And Embed key)
    HashFromPosition<<<grid,block>>>(d_p,wbl.d_HashKey,nElem,gC,hL);
    cudaDeviceSynchronize();
    
    CHECK(cudaGetLastError());



    
    //Kernel:Sort Key based on Hash.
    int nLog2Elem = (int)log2f((float)(nElem)+0.1f);
    //printf ("%d",nLog2Elem);
    sort(wbl.d_HashKey,wbl.d_HashKeyWork,nLog2Elem); //TYUI! secound parameter is for only work buffer(x:in-out,o:in&out-work)
 
    
    CHECK(cudaGetLastError());
    
    //Kernel:Generate Hash range reference map(HRRM).
    cudaMemset(wbl.d_HRRM,0,nHash*2*sizeof(int));
    GenerateHRRM<<<grid,block>>>(wbl.d_HashKey,wbl.d_HRRM,nElem);
    cudaDeviceSynchronize();
    
    CHECK(cudaGetLastError());
    
    cudaMemcpy(wbl.h_HRRMTest,wbl.d_HRRM,nHash*2*sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(wbl.h_HashKeyTest,wbl.d_HashKey,nElem*2*sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(wbl.h_nElemF,d_p,nBytes*3,cudaMemcpyDeviceToHost);
    
    {
        int* ct;
        int* h_h = wbl.h_HashKeyTest;
        float* h_n = wbl.h_nElemF;
        ct = (int * )malloc(nElem*sizeof(int));
        memset(ct,0,nElem*sizeof(int));
        for (int j=0;j<nElem;j++){
            int h = h_h[j];
            int k = h_h[j+nElem];
            ct[k]++;
        }
        for (int j=0;j<nElem;j++){
            if (ct[j]!=1){
                printf("(%d,%d)",j,ct[j]);
            }
        }

        
    }


    
    //Kernel:using Non-Aligned data and HRRM and Sorted Key-Hash,calculate Force fast.
    CalculateForce_GPUSort<<<grid,block>>>(d_f,d_p,wbl.d_HRRM,wbl.d_HashKey,nElem,length,gC,hL);
    //CalculateForce_GPUNaive<<<grid,block>>>(d_f,d_p,nElem,length);
    
    CHECK(cudaGetLastError());
    
    cudaMemcpy(h_f,d_f,nBytes*3,cudaMemcpyDeviceToHost);
    
    CHECK(cudaGetLastError());

}


void cuMain(void (*grpc)(V3Buf buf) ){
    initRand();
    
    //Buffer Initialization (CPU)

    int nElem = 256*8*8;
    int nBytes = nElem * sizeof(float);
    float *h_p,*h_v,*h_f,*h_fd;
    h_p = (float*)malloc(nBytes*3);
    h_v = (float*)malloc(nBytes*3);
    h_f = (float*)malloc(nBytes*3);
    h_fd= (float*)malloc(nBytes*3);

    //Buffer Initialization (GPU)
    float *d_p,*d_f;
    int *d_HashKeyIn,*d_HashKeyOut,*d_HRRM;
    cudaMalloc(&d_p,nBytes*3);
    cudaMalloc(&d_f,nBytes*3);
    cudaMalloc(&d_HashKeyIn,nElem*2*sizeof(int));
    cudaMalloc(&d_HashKeyOut,nElem*2*sizeof(int));
    cudaMalloc(&d_HRRM,(1<<(HSF*3))*2*sizeof(int));//HRRM size is determined by HashMax * 2
    float *h_df; //Test Buf;
    int *h_dh,*h_d_HRRM; // Test Buf;
    h_df = (float *)malloc(nBytes*3);
    h_dh = (int *)malloc(nElem*2*sizeof(int));
    h_d_HRRM = (int *)malloc((1<<(HSF*3))*2*sizeof(int));

    WorkBufList wbl = {d_HashKeyIn,d_HashKeyOut,d_HRRM,h_dh,h_d_HRRM,h_df};

    //Buffer Setting

    float length;
    V3Buf h_v3pos = CreateUniformParticles(h_p,0.5f,nElem,&length);
    V3Buf h_v3vel = CreateRandomVelocity(h_v,nElem);
    
    printf("length%f\n",length);
    
    for (int i=0;i<nElem*3;i++){
        //h_v[i]*=10.0f;
        h_v[i]=(float)((i*7)%13)/13.0f*2.0f-1.0f;
        h_v[i]*=2.0f;
    }

    double potential;
    CalculateForce_UseGPU_Naive(h_p,h_f,d_p,d_f,nElem,length);

    float dt = 0.005;
    int it = 0;
    while(true){
    
        //Graphics Functon(External) :transfer postion buffer to graphics function;
        if (it%20==0) (*grpc)(h_v3pos);
        
        //Position Update
        for (int i=0;i<nElem*3;i++){
            float p = h_p[i];
            p+=dt*(h_v[i]+0.5f*dt*h_f[i]);
            p = p- floorf(p/length)*length;
            h_p[i] = p;
        }
        
        //Force buffer becomes old because of position updated.
        SwapFloatPointer(&h_f,&h_fd);
        
        //CalculateForce_UseGPU_Naive(h_p,h_f,d_p,d_f,nElem,length);
        
        CalculateForce_UseGPU(h_p,h_f,d_p,d_f,nElem,length,wbl);
        
        for (int i=0;i<nElem*3;i++){
            h_v[i]+=dt*0.5*(h_f[i]+h_fd[i]);
        }
        
        printf("%f",h_p[1000]);
        
        //CalculateHamiltonian(h_p,h_v,nElem,potential);
        
        it++;
    }
    
}