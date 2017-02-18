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

__host__ __device__ inline int getHash(float x,float y,float z,int gC,float hL){
    // gC is grid Count(1dim) 
    // hL is grid Length(one grid,1dim)
    int Hx = (int)floorf(x/hL); //!!! if x == hL * gC ,it cause error!!!
    int Hy = (int)floorf(y/hL);
    int Hz = (int)floorf(z/hL);
    return Hx + Hy*(1<<6) + Hz*(1<<12); //this implemention uses gC < (1<<6 or 64);
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

typedef struct{
    int *d_HashKeyIn; //nElem * 2 (pair of hash and key);
    int *d_HashKeyOut;
    int *d_HRRM;//?
    int *h_HashKeyTest;//TestBuf
} WorkBufList;


void CalculateForce_UseGPU(float* h_p,float *h_f,float* d_p,float *d_f,int nElem,float length,WorkBufList wbl){
    float cutRadius = 4.0f;
    int gC = floor(length/cutRadius); //grid Size(1dim)
    float hL = length/gC; //cutRadius
    
    
    //Need Memcpy h2d;
    //N I
    
    dim3 block(256);
    dim3 grid(nElem/block.x);

    //Kernel:Hashing (Generate Hash And Embed key)
    HashFromPosition<<<grid,block>>>(d_p,wbl.d_HashKeyIn,nElem,gC,hL);
    cudaDeviceSynchronize();
    
    cudaMemcpy(wbl.h_HashKeyTest,wbl.d_HashKeyIn,nElem*2*sizeof(int),cudaMemcpyDeviceToHost);
    

    int cor = 0;
    for (int i=0;i<nElem;i++){
        float x,y,z;
        int j = i;
        x = h_p[j];
        y = h_p[j+nElem];
        z = h_p[j+2*nElem];
        bool check = wbl.h_HashKeyTest[i]==getHash(x,y,z,gC,hL);
        if (!check)printf("E:%f,%f,%f,%d,%d,%d,%d,%d\n",x,y,z,getHash(x,y,z,gC,hL),wbl.h_HashKeyTest[i],(int)floor(x/hL),(int)floor(y/hL),(int)floor(z/hL));
        if (check)cor++;
    }
    printf("%d\n",cor);
    
    //Kernel:Sort Key based on Hash.
    //Kernel:Generate Hash range reference map(HRRM).
    //Kernel:using Non-Aligned data and HRRM and Sorted Key-Hash,calculate Force fast.
}


void cuMain(void (*grpc)(V3Buf buf) ){
    initRand();
    
    //Buffer Initialization (CPU)

    int nElem = 256*8;
    int nBytes = nElem * sizeof(float);
    float *h_p,*h_v,*h_f,*h_fd;
    h_p = (float*)malloc(nBytes*3);
    h_v = (float*)malloc(nBytes*3);
    h_f = (float*)malloc(nBytes*3);
    h_fd= (float*)malloc(nBytes*3);

    //Buffer Initialization (GPU)
    float *d_p,*d_f;
    int *d_HashKeyIn,*d_HashKeyOut;
    cudaMalloc(&d_p,nBytes*3);
    cudaMalloc(&d_f,nBytes*3);
    cudaMalloc(&d_HashKeyIn,nElem*2*sizeof(int));
    cudaMalloc(&d_HashKeyOut,nElem*2*sizeof(int));
    float *h_df; //Test Buf;
    int *h_dh; // Test Buf;
    h_df = (float *)malloc(nBytes*3);
    h_dh = (int *)malloc(nElem*2*sizeof(int));

    WorkBufList wbl = {d_HashKeyIn,d_HashKeyOut,nullptr,h_dh};

    //Buffer Setting

    float length;
    V3Buf h_v3pos = CreateUniformParticles(h_p,1.0f,nElem,&length);
    V3Buf h_v3vel = CreateRandomVelocity(h_v,nElem);
    
    printf("length%f",length);
    
    for (int i=0;i<nElem;i++){
        h_v[i]*=10.0f;
    }

    double potential;
    CalculateForce_UseGPU_Naive(h_p,h_f,d_p,d_f,nElem,length);

    float dt = 0.005;
    int it = 0;
    while(true){
    
        //Graphics Functon(External) :transfer postion buffer to graphics function;
        if (it%1==0) (*grpc)(h_v3pos);
        
        //Position Update
        for (int i=0;i<nElem*3;i++){
            float p = h_p[i];
            p+=dt*(h_v[i]+0.5f*dt*h_f[i]);
            p = p- floorf(p/length)*length;
            h_p[i] = p;
        }
        
        //Force buffer becomes old because of position updated.
        SwapFloatPointer(&h_f,&h_fd);
        
        CalculateForce_UseGPU_Naive(h_p,h_f,d_p,d_f,nElem,length);
        
        CalculateForce_UseGPU(h_p,h_f,d_p,d_f,nElem,length,wbl);
        
        for (int i=0;i<nElem*3;i++){
            h_v[i]+=dt*0.5*(h_f[i]+h_fd[i]);
        }
        
        //CalculateHamiltonian(h_p,h_v,nElem,potential);
        
        it++;
    }
    
}