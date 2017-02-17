#include <stdio.h>
#include <time.h>
#include "cuHeader.h"
#include <math.h>

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
V3Buf CreateUniformParticles(float* buf,float rho,int nElem){
    int nloop = (int)powf((float)(nElem/4)+0.1f,1.0f/3);
    printf("%d\n",nloop);
    
    float volume = powf(nElem/rho,1.0f/3);
    float a = volume/nloop;
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

V3Buf CalculateForce(float *force,float *pos,int nElem){
    float *h_fx,*h_fy,*h_fz,*h_px,*h_py,*h_pz;
    h_fx = force;
    h_fy = h_fx + nElem;
    h_fz = h_fy + nElem;
    h_px = pos;
    h_py = h_px + nElem;
    h_pz = h_py + nElem;
    
    for(int i=0;i<nElem;i++){
        float px=h_px[i],py=h_py[i],pz=h_pz[i];
        float r2 = px*px+py*py+pz*pz;
        float fx,fy,fz,k;
        k = 10.0f;
        fx = -px/r2*k;
        fy = -py/r2*k;
        fz = -pz/r2*k;
        h_fx[i] = fx;
        h_fy[i] = fy;
        h_fz[i] = fz;
    }
    V3Buf v3Buf = {h_fx,h_fy,h_fz,nElem};
    return v3Buf;
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
    
    V3Buf h_pos = CreateUniformParticles(h_p,0.5f,nElem);
    V3Buf h_vel = CreateRandomVelocity(h_v,nElem);
    
    for (int i=0;i<nElem;i++){
            h_v[i]+=4.0;
    }

    CalculateForce(h_f,h_p,nElem);//Zeroth Calculation;

    float dt = 0.01;
    int it = 0;
    while(true){
    
        //Graphics Functon(External) :transfer postion buffer to graphics function;
        if (it%30==0) (*grpc)(h_pos);
        
        for (int i=0;i<nElem*3;i++){
            h_p[i]+=dt*(h_v[i]+0.5f*dt*h_f[i]);
        }
        
        SwapFloatPointer(&h_f,&h_fd);
        CalculateForce(h_f,h_p,nElem);
        
        for (int i=0;i<nElem*3;i++){
            h_v[i]+=dt*0.5*(h_f[i]+h_fd[i]);
        }
        
        it++;
    }
    
}