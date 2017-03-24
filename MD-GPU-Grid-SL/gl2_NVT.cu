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

typedef struct{
    float Lx;
    float Ly;
    float Lz;
    int gCx;
    int gCy;
    int gCz;
    float hLx;
    float hLy;
    float hLz;
    float dt;
    int nElem;
    double sumPotential;
    double sumKinetic;
    double sumForceDotPos;
}SystemInfo;

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
    for(int i=0;i<nElem;i++){
        float px=h_px[i],py=h_py[i],pz=h_pz[i];
        float vx=h_vx[i],vy=h_vy[i],vz=h_vz[i];
        float r = sqrtf(px*px+py*py+pz*pz);
        float v2= vx*vx+vy*vy+vz*vz;
        energy += (double)(v2/2);
    }
    printf("%lf  %lf  %lf\n",energy+potential,energy,potential);
    
    return energy+potential;
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
    //double t_fx=0.0f,t_fy=0.0f,t_fz=0.0f;
    
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


__global__ void UpdateDifferentialMomentum(float *f,float *v,int nElem,float dt,float s){
 
    int idx = threadIdx.x + blockIdx.x * blockDim.x; //1dim grid ,1im block
    
    float tmp = v[idx]*s;
    tmp += dt*0.5f*f[idx];
    v[idx] = tmp;
 
 
 /*
         for (int i=0;i<nElem*3;i++){
            h_v[i]+=dt*0.5*h_f[i];
        }
        */   
}

__global__ void UpdateDifferentialPosition(float *p,float *v,int nElem,float dt,SystemInfo *sysInfo){

    int idx = threadIdx.x + blockIdx.x * blockDim.x; //1dim grid ,1im block

    float length;
    length = idx/nElem == 0 ? sysInfo->Lx : length;
    length = idx/nElem == 1 ? sysInfo->Ly : length;
    length = idx/nElem == 2 ? sysInfo->Lz : length;

    float tmp = p[idx];
    tmp+=dt*v[idx];
    tmp = tmp- floorf(tmp/length)*length;
    p[idx] = tmp;

/*
            float p = h_p[i];
            p+=dt*h_v[i];
            p = p- floorf(p/length)*length;
            h_p[i] = p;*/
}

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
    __syncthreads();
    
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


__device__ double __atomicAdd(double* address, double val)
{
  unsigned long long int* address_as_ull =
    (unsigned long long int*)address;
  unsigned long long int old = *address_as_ull, assumed;
  do {
    assumed = old;
    old = atomicCAS(address_as_ull, assumed,
		    __double_as_longlong(val +
					 __longlong_as_double(assumed)));
  } while (assumed != old);
  return __longlong_as_double(old);
}

__global__ void CalculateSumOfPotential(float *potential,SystemInfo *sysInfo){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int sdx = threadIdx.x;
    //int wdx = threadIdx.x % 32;
    
    extern __shared__ float smem[];
    smem[sdx] = potential[idx];
    __syncthreads();
    
    double sum = (double)blockLevelSum(smem);
    if (sdx==0)__atomicAdd(&(sysInfo->sumPotential),sum/2);
    //２体間ポテンシャルをダブルでカウントしてしまっているため２で割る！！
}

__global__ void CalculateSumOfForceDotPos(float *fdp,SystemInfo *sysInfo){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int sdx = threadIdx.x;
    //int wdx = threadIdx.x % 32;
    
    extern __shared__ float smem[];
    smem[sdx] = fdp[idx];
    __syncthreads();
    
    double sum = (double)blockLevelSum(smem);
    if (sdx==0)__atomicAdd(&(sysInfo->sumForceDotPos),sum/2);
    //２体間計算をダブルでカウントしてしまっているため２で割る！！
}

__global__ void CalculateSumOfKinetic(float *velocity,SystemInfo *sysInfo){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int sdx = threadIdx.x;
    //int wdx = threadIdx.x % 32;
    
    extern __shared__ float smem[];
    float tmp = velocity[idx];
    smem[sdx] = tmp * tmp;
    __syncthreads();
    
    double sum = (double)blockLevelSum(smem);
    if (sdx==0)__atomicAdd(&(sysInfo->sumKinetic),sum/2);
    //運動エネルギー計算の係数が1/2なので２で割る！！
}

#define HSF 6//hash size = 1<<HSF

__host__ __device__ inline int getHashPartX(float v,SystemInfo sysInfo){
    float hLx = sysInfo.hLx;
    return (int)floorf(v/hLx);
}

__host__ __device__ inline int getHashPartY(float v,SystemInfo sysInfo){
    float hLy = sysInfo.hLy;
    return (int)floorf(v/hLy);
}

__host__ __device__ inline int getHashPartZ(float v,SystemInfo sysInfo){
    float hLz = sysInfo.hLz;
    return (int)floorf(v/hLz);
}

__host__ __device__ inline int HashParts2Hash(int Hx,int Hy,int Hz){
    return Hx + Hy*(1<<HSF) + Hz*(1<<(2*HSF)); //this implemention uses gC < (1<<HSF);
}

__host__ __device__ inline int getHash(float x,float y,float z,SystemInfo sysInfo){
    int Hx = getHashPartX(x,sysInfo); //!!! if x == hL * gC ,it cause error!!!
    int Hy = getHashPartY(y,sysInfo);
    int Hz = getHashPartZ(z,sysInfo);
    return HashParts2Hash(Hx,Hy,Hz);
}

__host__ __device__ inline void recoverHash(int hash,int &x,int &y,int &z){
    //Not Implement;
}

__global__ void HashFromPosition(float *d_p,int *d_hash,SystemInfo *sysInfo){
    int idx = threadIdx.x + blockIdx.x * blockDim.x;
    int nElem = sysInfo->nElem;
    float x = d_p[idx];
    float y = d_p[nElem+idx];
    float z = d_p[2*nElem+idx];
    int h = getHash(x,y,z,*sysInfo);
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

__global__ void AlignMemory(float *t,float *s,int *key,int nElem){
    int idx = threadIdx.x + blockIdx.x * blockDim.x; //1dim grid ,1im block
    t[idx] = s[key[idx]];
}

__device__ inline void warpLevelIndex(int* s,int *t){
    int sd = threadIdx.x;
    int wd = sd%32;
    
    int tmp,swp;
    tmp = *s;
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
    
    *t = tmp;
}

__global__ void CalculateForce_GPUSort(float *force,float *pos,int* hrrm,SystemInfo *sysInfo,float *d_pot,float *d_fdp){
    int idx = threadIdx.x + blockIdx.x * blockDim.x; //1dim grid ,1im block
    int sdx = threadIdx.x;
    int wdx = sdx%32;
    int nElem = sysInfo->nElem;
    float Lx = sysInfo->Lx;
    float Ly = sysInfo->Ly;
    float Lz = sysInfo->Lz;
    int gCx = sysInfo->gCx;
    int gCy = sysInfo->gCy;
    int gCz = sysInfo->gCz;
    float *d_fx,*d_fy,*d_fz,*d_px,*d_py,*d_pz;
    d_fx = force;
    d_fy = d_fx + nElem;
    d_fz = d_fy + nElem;
    d_px = pos;
    d_py = d_px + nElem;
    d_pz = d_py + nElem;

        
    extern __shared__ float smem[];
    
    
    //SYOKIKA
    float t_fx=0.0f,t_fy=0.0f,t_fz=0.0f;
    //double t_fx=0.0f,t_fy=0.0f,t_fz=0.0f;
    
    double potential = 0.0;// must implement to calculate Hamiltoniam...
    double forceDotPos = 0.0;// must implement to calculate Pressure...
    
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
    int Hx = getHashPartX(px,*sysInfo);
    int Hy = getHashPartY(py,*sysInfo);
    int Hz = getHashPartZ(pz,*sysInfo);
    int Hh = HashParts2Hash(Hx,Hy,Hz);
    int nHash = 1<<(HSF*3);
    
    int hid;
    warpLevelIndex(&Hh,&hid);
    
    //if (blockIdx.x == 30 & threadIdx.x/32 == 0)printf("<%d,%d>",threadIdx.x,hid);
    
    bool warpSync = __all(hid == 0);
    //if (threadIdx.x%32 == 31 & warpSync )printf("[!]");
    
    
    for (int c=0;c<27;c++){
        int dHx = c%3-1;
        int dHy = (c/3)%3-1;
        int dHz = (c/3/3)%3-1;
        int pHx = (Hx+gCx+dHx)%gCx;
        int pHy = (Hy+gCy+dHy)%gCy;
        int pHz = (Hz+gCz+dHz)%gCz;
        int h = HashParts2Hash(pHx,pHy,pHz);
        int start = hrrm[h];
        int end = hrrm[h+nHash];
        if(!(start<=end & end-start<1000000)){printf("D");};
        
        for (int k=start;k<end;k++){
            //int j = hash[k+nElem];
            //if (hash[k]!=h)printf("H");
            int j = (k)%(end-start) + start;
            if (i==j)continue;
            float dx,dy,dz,r2,r2i,r06i,r12i,fc,fx,fy,fz;
            float qx,qy,qz;
            qx = d_px[j];
            qy = d_py[j];
            qz = d_pz[j];
            dx = px-qx;
            dy = py-qy;
            dz = pz-qz;
            if(dx<-Lx/2) dx+=Lx;
            if(dx> Lx/2) dx-=Lx;
            if(dy<-Ly/2) dy+=Ly;
            if(dy> Ly/2) dy-=Ly;
            if(dz<-Lz/2) dz+=Lz;
            if(dz> Lz/2) dz-=Lz;
            //if (c==13 & dx>hLx*3)printf("G");
            if (dx > Lx*2 ) printf("P");
            if (dy > Ly*2 ) printf("P");
            if (dz > Lz*2 ) printf("P");
            if (getHash(qx,qy,qz,*sysInfo)!=h)printf("E");
            r2 = dx*dx+dy*dy+dz*dz;
            r2i = 1.0f/r2;
            r06i = r2i * r2i * r2i;
            r12i = r06i * r06i;
            potential += (ce12 * r12i - ce06 * r06i);
            fc = (cf12*r12i-cf06*r06i)*r2i;
            fx = fc * dx;
            fy = fc * dy;
            fz = fc * dz;
            forceDotPos -= r2 * fc;
            t_fx+=fx;
            t_fy+=fy;
            t_fz+=fz;
        }        
        
    }
        
    d_fx[i]=t_fx;
    d_fy[i]=t_fy;
    d_fz[i]=t_fz;
    d_pot[i]=potential;
    d_fdp[i]= forceDotPos;
}

//規約　引数渡しとしては使わない．only for「WorkBuffer」
typedef struct{ //CPU level struct
    float * h_pot;
    int *d_HashKey; //nElem * 2 (pair of hash and key);
    int *d_HashKeyWork;
    int *d_HRRM;//?
    int *h_HashKeyTest;//TestBuf
    int *h_HRRMTest;
    float *h_nElemF;
    float *h_f;
    float *h_p;
    float *h_v;
    SystemInfo *h_sysInfo;
} WorkBufList;


void CalculateForce_UseGPU(float *&d_f,float* &d_p,float * &d_v,float* &d_p2,float * &d_v2,float* &d_pot,float* &d_fdp,SystemInfo *d_sysInfo,WorkBufList wbl){
    cudaMemcpy(wbl.h_sysInfo,d_sysInfo,sizeof(SystemInfo),cudaMemcpyDeviceToHost);
    
    int nElem = wbl.h_sysInfo->nElem;
    float Lx = wbl.h_sysInfo->Lx;
    float Ly = wbl.h_sysInfo->Ly;
    float Lz = wbl.h_sysInfo->Lz;
    float cutRadius = 4.0f;
    int gCx = floor(Lx/cutRadius); //grid Size(1dim)
    int gCy = floor(Ly/cutRadius);
    int gCz = floor(Lz/cutRadius);
    float hLx = Lx/gCx; //cutRadius
    float hLy = Ly/gCy;
    float hLz = Lz/gCz;
    int nHash = 1<<(HSF*3);
    int nBytes = nElem * sizeof(float);
    
    wbl.h_sysInfo-> gCx = gCx;
    wbl.h_sysInfo-> gCy = gCy;
    wbl.h_sysInfo-> gCz = gCz;
    wbl.h_sysInfo-> hLx = hLx;
    wbl.h_sysInfo-> hLy = hLy;
    wbl.h_sysInfo-> hLz = hLz;
    cudaMemcpy(d_sysInfo,wbl.h_sysInfo,sizeof(SystemInfo),cudaMemcpyHostToDevice);
    
    
    dim3 block(128);
    dim3 grid(nElem/block.x);
    //printf("G<%d>",grid.x);
    
    CHECK(cudaGetLastError());

    //Kernel:Hashing (Generate Hash And Embed key)
    HashFromPosition<<<grid,block>>>(d_p,wbl.d_HashKey,d_sysInfo);
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
    

    AlignMemory<<<grid,block>>>(d_p2+nElem*0,d_p+nElem*0,wbl.d_HashKey+nElem,nElem);
    AlignMemory<<<grid,block>>>(d_p2+nElem*1,d_p+nElem*1,wbl.d_HashKey+nElem,nElem);
    AlignMemory<<<grid,block>>>(d_p2+nElem*2,d_p+nElem*2,wbl.d_HashKey+nElem,nElem);
    AlignMemory<<<grid,block>>>(d_v2+nElem*0,d_v+nElem*0,wbl.d_HashKey+nElem,nElem);
    AlignMemory<<<grid,block>>>(d_v2+nElem*1,d_v+nElem*1,wbl.d_HashKey+nElem,nElem);
    AlignMemory<<<grid,block>>>(d_v2+nElem*2,d_v+nElem*2,wbl.d_HashKey+nElem,nElem);
    SwapFloatPointer(&d_p,&d_p2);
    SwapFloatPointer(&d_v,&d_v2);
    cudaDeviceSynchronize();
    
    CHECK(cudaGetLastError());
    
    //Kernel:using Non-Aligned data and HRRM and Sorted Key-Hash,calculate Force fast.
    CalculateForce_GPUSort<<<grid,block,block.x*sizeof(float)*3>>>(d_f,d_p,wbl.d_HRRM,d_sysInfo,d_pot,d_fdp);
    //CalculateForce_GPUNaive<<<grid,block>>>(d_f,d_p,nElem,length);
    
   
    
    CHECK(cudaGetLastError());
    
    /*
    cudaMemcpy(wbl.h_HRRMTest,wbl.d_HRRM,nHash*2*sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(wbl.h_HashKeyTest,wbl.d_HashKey,nElem*2*sizeof(int),cudaMemcpyDeviceToHost);
    cudaMemcpy(wbl.h_nElemF,d_p,nBytes*3,cudaMemcpyDeviceToHost);
    cudaMemcpy(wbl.h_f,d_f,nBytes*3,cudaMemcpyDeviceToHost);
    
    cudaMemcpy(wbl.h_pot,wbl.d_pot,nBytes,cudaMemcpyDeviceToHost);
    cudaMemcpy(wbl.h_p,d_p,nBytes*3,cudaMemcpyDeviceToHost);
    cudaMemcpy(wbl.h_v,d_v,nBytes*3,cudaMemcpyDeviceToHost);
    */
    CHECK(cudaGetLastError());

}


void cuMain(void (*grpc)(V3Buf buf) ){
    initRand();
    
    //Buffer Initialization (CPU)

    int nElem = 256*8*8;
    int nBytes = nElem * sizeof(float);
    float *h_p,*h_v,*h_f,*h_pot;
    SystemInfo *h_sysInfo;
    h_p = (float*)malloc(nBytes*3);
    h_v = (float*)malloc(nBytes*3);
    h_f = (float*)malloc(nBytes*3);
    h_pot = (float*)malloc(nBytes);
    h_sysInfo = (SystemInfo*)malloc(sizeof(SystemInfo));

    //Buffer Initialization (GPU)
    float *d_p,*d_f,*d_pot,*d_v,*d_p2,*d_v2,*d_fdp;
    int *d_HashKeyIn,*d_HashKeyOut,*d_HRRM;
    SystemInfo *d_sysInfo;
    cudaMalloc(&d_p,nBytes*3);
    cudaMalloc(&d_f,nBytes*3);
    cudaMalloc(&d_v,nBytes*3);
    cudaMalloc(&d_HashKeyIn,nElem*2*sizeof(int));
    cudaMalloc(&d_HashKeyOut,nElem*2*sizeof(int));
    cudaMalloc(&d_HRRM,(1<<(HSF*3))*2*sizeof(int));//HRRM size is determined by HashMax * 2
    cudaMalloc(&d_pot,nBytes);
    cudaMalloc(&d_p2,nBytes*3);
    cudaMalloc(&d_v2,nBytes*3);
    cudaMalloc(&d_fdp,nBytes);
    cudaMalloc((void**)&d_sysInfo,sizeof(SystemInfo));
    float *h_df; //Test Buf;
    int *h_dh,*h_d_HRRM; // Test Buf;
    h_df = (float *)malloc(nBytes*3);
    h_dh = (int *)malloc(nElem*2*sizeof(int));
    h_d_HRRM = (int *)malloc((1<<(HSF*3))*2*sizeof(int));

    WorkBufList wbl = 
        {h_pot,d_HashKeyIn,d_HashKeyOut,d_HRRM,h_dh,h_d_HRRM,h_df,h_f,h_p,h_v,h_sysInfo};

    //Buffer Setting

    float length;
    //V3Buf h_v3pos = CreateUniformParticles(h_p,1.0f,nElem,&length);
    V3Buf h_v3pos = CreateUniformParticles(h_p,0.42f,nElem,&length);
    V3Buf h_v3vel = CreateRandomVelocity(h_v,nElem);
    
    length *= 1;
    
    printf("length%f\n",length);
    
    for (int i=0;i<nElem*3;i++){
        //h_v[i]*=10.0f;
        h_v[i]=(float)((i*7)%13)/13.0f*2.0f-1.0f;
        h_v[i]*=4.0f;//2.3f;
    }


    float dt = 0.01 * 0.25;

    
    h_sysInfo->Lx = length * 1.0f;//*1.1f -> L-S boundary???
    h_sysInfo->Ly = length * 1.0f;
    h_sysInfo->Lz = length * 1.0f;
    h_sysInfo->dt = dt;
    h_sysInfo->nElem = nElem;
    cudaMemcpy(d_sysInfo,h_sysInfo,sizeof(SystemInfo),cudaMemcpyHostToDevice);    
    
    
    cudaMemcpy(d_p,h_p,nBytes*3,cudaMemcpyHostToDevice);
    cudaMemcpy(d_v,h_v,nBytes*3,cudaMemcpyHostToDevice);
    
    double T_target = 2.0;
    double s = 1.0;
    
    ///// first calc : F and Pot
    CalculateForce_UseGPU(d_f,d_p,d_v,d_p2,d_v2,d_pot,d_fdp,d_sysInfo,wbl);
    
    h_sysInfo->sumPotential = 0.0;
    cudaMemcpy(d_sysInfo,h_sysInfo,sizeof(SystemInfo),cudaMemcpyHostToDevice);
    CalculateSumOfPotential<<<nElem/1024,1024,1024*sizeof(float)>>>(d_pot,d_sysInfo);
    cudaDeviceSynchronize();
    cudaMemcpy(h_sysInfo,d_sysInfo,sizeof(SystemInfo),cudaMemcpyDeviceToHost);
    ////

    int it = 0;
    while(true & it * dt < 3000){
    
        
        dim3 block(256);
        dim3 grid(nElem*3/block.x);
       
        UpdateDifferentialMomentum<<<grid,block>>>(d_f,d_v,nElem,dt,1);
        cudaDeviceSynchronize();
        
        //Calculation KineticEnergy
        h_sysInfo->sumKinetic = 0.0;
        cudaMemcpy(d_sysInfo,h_sysInfo,sizeof(SystemInfo),cudaMemcpyHostToDevice);
        CalculateSumOfKinetic<<<nElem*3/1024,1024,1024*sizeof(float)>>>(d_v,d_sysInfo);
        cudaDeviceSynchronize();
        cudaMemcpy(h_sysInfo,d_sysInfo,sizeof(SystemInfo),cudaMemcpyDeviceToHost);
        
        UpdateDifferentialPosition<<<grid,block>>>(d_p,d_v,nElem,dt,d_sysInfo);
        cudaDeviceSynchronize();
        
        CalculateForce_UseGPU(d_f,d_p,d_v,d_p2,d_v2,d_pot,d_fdp,d_sysInfo,wbl);
        cudaDeviceSynchronize();
        
        h_sysInfo->sumPotential = 0.0;
        cudaMemcpy(d_sysInfo,h_sysInfo,sizeof(SystemInfo),cudaMemcpyHostToDevice);
        CalculateSumOfPotential<<<nElem/1024,1024,1024*sizeof(float)>>>(d_pot,d_sysInfo);
        cudaDeviceSynchronize();
        cudaMemcpy(h_sysInfo,d_sysInfo,sizeof(SystemInfo),cudaMemcpyDeviceToHost);
        
        h_sysInfo->sumForceDotPos = 0.0;
        cudaMemcpy(d_sysInfo,h_sysInfo,sizeof(SystemInfo),cudaMemcpyHostToDevice);
        CalculateSumOfForceDotPos<<<nElem/1024,1024,1024*sizeof(float)>>>(d_fdp,d_sysInfo);
        cudaDeviceSynchronize();
        cudaMemcpy(h_sysInfo,d_sysInfo,sizeof(SystemInfo),cudaMemcpyDeviceToHost);       
        
        {
            double T = h_sysInfo->sumKinetic/nElem/3 * 2;
            s = sqrt(T_target/T);
            if (it * dt > 5) s = 1.0;
        }
        
        
        UpdateDifferentialMomentum<<<grid,block>>>(d_f,d_v,nElem,dt,s);
        cudaDeviceSynchronize();
        
        //Graphics Functon(External) :transfer postion buffer to graphics function;
        if (it%100==1) {;//20ステップに１回表示：粒子数が多いと表示で律速するので．
           
            cudaMemcpy(wbl.h_p,d_p,nBytes*3,cudaMemcpyDeviceToHost);
            cudaMemcpy(wbl.h_v,d_v,nBytes*3,cudaMemcpyDeviceToHost);            
            
            (*grpc)(h_v3pos);
            
            double V =  h_sysInfo->Lx * h_sysInfo->Ly * h_sysInfo->Lz;
            double rho = ( (double)nElem ) / V;
            double H = h_sysInfo->sumPotential+h_sysInfo->sumKinetic;
            double T = h_sysInfo->sumKinetic/nElem/3 * 2;
            double P = ( T * nElem - h_sysInfo->sumForceDotPos/3.0 )/V;
            printf("%lf,T=%lf,P=%lf,rho=%lf:gpu\n",H,T,P,rho);
            printf("t:%f\n",(float)(it)*dt);
        }
        
        
        
        it++;
    }
    
}