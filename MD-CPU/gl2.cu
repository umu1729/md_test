#include <stdio.h>
#include <time.h>
#include "cuHeader.h"
#include <math.h>
#include <cuda_runtime.h>

void initRand(){
    time_t t;
    srand((unsigned) time(&t));

}

void initData(float *ip,int size){
    for (int i=0;i<size;i++){
        ip[i] = (float)( rand() & 0xFF )/256*10;
    }
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
    float k = 10.0f;
    
    double energy = 0.0;
    for(int i=0;i<1;i++){
        float px=h_px[i],py=h_py[i],pz=h_pz[i];
        float vx=h_vx[i],vy=h_vy[i],vz=h_vz[i];
        float r = sqrtf(px*px+py*py+pz*pz);
        float v2= vx*vx+vy*vy+vz*vz;
        energy += (double)(v2/2);
    }
    energy += potential;
    printf("%lf\n",energy);
    
    return energy;
}

void SwapFloatPointer(float** a,float **b){
    float *tmp = *b;
    *b = *a;
    *a = tmp;
}

void cuMain(void (*grpc)(V3Buf buf) ){
    initRand();

    int nElem = 256*8;
    int nBytes = nElem * sizeof(float);
    float *h_p,*h_v,*h_f,*h_fd;
    h_p = (float*)malloc(nBytes*3);
    h_v = (float*)malloc(nBytes*3);
    h_f = (float*)malloc(nBytes*3);
    h_fd= (float*)malloc(nBytes*3);
    //cudaMallocHost(&h_fd,nBytes*3);
    
    float length;
    V3Buf h_v3pos = CreateUniformParticles(h_p,1.5f,nElem,&length);
    V3Buf h_v3vel = CreateRandomVelocity(h_v,nElem);
    
    for (int i=0;i<nElem;i++){
        h_v[i]*=10.0f;
    }

    double potential;
    CalculateForce(h_f,h_p,nElem,length,&potential);//Zeroth Calculation;

    float dt = 0.005;
    int it = 0;
    while(true){
    
        //Graphics Functon(External) :transfer postion buffer to graphics function;
        if (it%1==0) (*grpc)(h_v3pos);
        
        for (int i=0;i<nElem*3;i++){
            float p = h_p[i];
            p+=dt*(h_v[i]+0.5f*dt*h_f[i]);
            p = p- floorf(p/length)*length;
            h_p[i] = p;
        }
        
        SwapFloatPointer(&h_f,&h_fd);
        
        
        CalculateForce(h_f,h_p,nElem,length,&potential);
        
        for (int i=0;i<nElem*3;i++){
            h_v[i]+=dt*0.5*(h_f[i]+h_fd[i]);
        }
        
        CalculateHamiltonian(h_p,h_v,nElem,potential);
        
        it++;
    }
    
}