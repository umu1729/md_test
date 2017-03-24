#pragma once
#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>
#include <string.h>

#include "../serializer/serial.h"


int main(int argc,char** argv){
    FILE* fp;
    fp = fopen("rawData/test.rawdata","rb");
    fseek(fp,0,SEEK_END);
    size_t size = ftell(fp);
    fseek(fp,0,SEEK_SET);
    printf("fsize:%d\n",(int)size);
    void *dat = malloc(size);
    fread(dat,1,size,fp);
    
    void *now = dat;
    void *tmp = now;
    float *target;
    int nElem;
    
    for (int j=0;j<5;j++){
        now = NextBuf(now);   
        tmp = DescendBuf(now);
        printf("%s\n",NameBuf(tmp));
        target = (float*)DescendBuf(tmp);
        nElem = SizeBuf(tmp)/sizeof(float);
        printf("<nELem:%d>",nElem);
        for (int i=0;i<10;i++){
            printf("%f,",target[i]);
        }
        printf("\n");
    }    
    
    return 0;
}