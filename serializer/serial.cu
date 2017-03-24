#pragma once
#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>
#include <string.h>

#include "serial.h"


int main(int argc,char** argv){


    //  �f�[�^��p�ӂ���

    int d1s = 128;
    int* d1 = (int*)malloc(d1s*sizeof(d1s));
    
    for (int i=0;i<d1s;i++){
        d1[i] = i;
    }
    
    int d2s = 64;
    float* d2 = (float*)malloc(d2s*sizeof(float));
    
    int d3s = 16;
    double* d3 = (double*) malloc(d3s*sizeof(double));
    
    for (int i=0;i<d3s;i++){
        d3[i] = ((double)i)/12;
    }
    
    //ChunkIndex�̍쐬
    //ChunkIndex�́C�ΏۂƂ���o�b�t�@�i�z��Ƃ��j�̃|�C���^��C����ChunkIndex�C�q��ChunkIndex��ێ����܂��D
    //ChunkIndex���쐬���邱�ƂŁC�f�[�^�\�����`���C���[�h���Ɏg���₷���o�͂��\�����܂��D
    // constuctChunk(name,size,pointer)�́Cpointer��擪�Ƃ��āC����size������ChunkIndex���C
    // name�Ƃ������O�ō쐬���܂��D
    // groupChunk(name,ChunkIndexPointerArray,count)�́C�O���[�v��������ChunkIndex�́C
    // �|�C���^�[�z��(ChunkIndexPointerArray)�ƁC���̔z��̒���count����Cname�Ƃ������O��
    // �O���[�v�����ꂽChunkIndex���쐬���܂�,����͐e�q�\�����\�z���܂��D
    
    ChunkIndex chk1 = constructChunk("int1",d1s*sizeof(int),d1);
    ChunkIndex chk2 = constructChunk("flo1",d2s*sizeof(float),d2);
    ChunkIndex* chksA[2] = {&chk1,&chk2};
    ChunkIndex chk3 = groupChunk("g1",chksA,2);
    ChunkIndex chk4 = constructChunk("dou1",d3s*sizeof(double),d3);
    ChunkIndex* chksB[2] = {&chk3,&chk4};
    ChunkIndex chk5 = groupChunk("g2",chksB,2);
    
    //�ȏ�̒�`�ŁCg2{g1{int1,flo1},dou1}�Ƃ����悤��ChunkIndex���쐬���ꂽ�D
    //�������C�L�@�@groupName{child1Name,child2Name,...}�ƂȂ��Ă���D
    // ChunkIndex���O���[�v�����邱�ƂŁC�o�̓t�@�C�����Ńf�[�^���K�w���i�e���q�����j�ł��܂��D
    
    // Chunk���̍X�V
    //updateChunkInfo(topLevelChunkIndex)�֐���p���āCChunkIndex�̏����ŐV�̏�Ԃɂ��܂��D
    //ChunkIndex��V���ɍ������C�O���[�v�������ꍇ�́C��ԏ��ChunkIndex(topLevelChunkIndex)��
    //�����ɉ����āD�����X�V����K�v������܂��D����͓����I�ɂ́C�t�@�C���o�͎��̃o�b�t�@�T�C�Y��
    //�v�Z���s���Ă��܂��D
    
    updateChunkInfo(chk5);
    
    //�o�͗p�t�@�C���n���h����p�ӂ��āCChunkIndex�����ɁCsaveChunk(fp,chuk)�֐��Ńf�[�^���������݂܂��D
    //ChunkIndex���O���[�v�����ꂽ���̂ł���΁C�O���[�v�S�̂��������܂�܂��D
    //��{�I�ɂ́C��ԏ�̐e�������ɓn���܂�
    //�i���ӁF���݂ł͈�ԏ�̐e���n����邱�Ƃ�z�肵�Ă���̂ŁC�q������Ǝv��ʃo�O�����邩������܂���D
    //fp�͗p�ӂ����t�@�C���n���h���Cchunk�͏�������ChunkIndex�ł��D
    //saveChunk���P�x�g�p���邲�ƂɃt�@�C���Ƀf�[�^���������܂�܂��D
    //���̃T���v���ł̓`�����N���Q�񏑂����܂�Ă��܂��D
    //���̏ꍇ�ł͂��܂�Ӗ�������܂��񂪁C�Ⴆ��saveChunk�̊Ԃ�ChunkIndex�ɂ���ĊĎ�����Ă���
    //�|�C���^�[���ύX���ꂽ�ꍇ�iMD�ł����1step�i�߂�Ȃǁj
    //saveChunk�ɂ���āC�قȂ�f�[�^��ۑ����邱�Ƃ��ł��܂��D
    //MD�̏ꍇ�́Cchunkindex�Ɉʒu���o�b�t�@�����Ă����C��step���Ƃ�saveChunk�����邱�ƂŁC
    //�t�@�C���Ɋe�X�e�b�v�̏�񂪂��ׂďo�͂ł��܂��D
    
    //�o�̓t�@�C���̍\���́C�P��ChunkIndex���ƂɁCHeader(�f�[�^�̏��)��Body(�f�[�^�{��)
    //���A�����ď������܂��\���ɂȂ��Ă��܂�(Header+Body��Chunk�ƌĂԂ��Ƃɂ��܂�)�D
    //�e�q���Ȃ��P���ChunkIndex��saveChunk�ŉ��x���������ޏꍇ�́C
    //Header-Body-Header-Body-Header-Body-...�Ƃ����\���ɂȂ��Ă��܂��D
    //�܂��́CChunk-Chunk-Chunk-...
    //�e�q�֌W�̂���ChunkIndex��ۑ������ꍇ�́C
    //�e��Body�����ɁC�q��Chunk�����ߍ��܂�܂��D
    //parent{child}�̏ꍇ�́C
    //Hp-Hc-Bc-Hp-Hc-Bc-... �ƂȂ�܂� Hp:parent's header,Bc:child's header
    //�q�������̏ꍇ�́C
    //Hp-Hc1-Bc1-Hc2-Bc2-Hp-Hc1-Bc1-Hc2-Bc2-..
    //�ƂȂ�܂��DHeader�ɁCBody�̃T�C�Y�̏�񂪏������܂�Ă���̂ŁC�ȒP�Ȋ֐����g����
    //�q��Header�⎟��Header�̃o�b�t�@�[�ʒu�����ł��܂��D
    
    FILE * fp= fopen("dat.txt","wb");
    printf("ch:%d\n",(int)sizeof(ChunkHeader));
    saveChunk(fp,chk5);
    saveChunk(fp,chk5);
    fclose(fp);
    
    //��������́C�o�͂��ꂽ�t�@�C������C�f�[�^�𕜌�����R�[�h�ł��D
    
    printf("LOAD\n");
    
    //�܂��ŏ��ɁC�o�̓t�@�C�������ׂēǂݍ��݂܂��D
    //���F���ׂēǂݍ��܂Ȃ��ŏ��X�ɓǂݍ��ޕ������J������\��ł��D
    
    fp = fopen("dat.txt","rb");
    fseek(fp,0,SEEK_END);
    size_t size = ftell(fp);
    fseek(fp,0,SEEK_SET);
    printf("fsize:%d\n",(int)size);
    void *dat = malloc(size);
    fread(dat,1,size,fp);
    
    //�ǂݍ��ނƁC�o�b�t�@�|�C���^���C��ԍŏ��ɏ������܂ꂽChunk�̐擪���w���Ă����ԂƂȂ�܂��D
    //(�܂�CHeader�̐擪�ł�����j
    //��`����Ă���֗��֐���p���āC�o�b�t�@�|�C���^���C�ړ������邱�Ƃɂ���ĎQ�Ƃł���悤�ɂ��܂��D
    //NameBuf(ChunkPointer)�́C�`�����N�̐擪�|�C���^����C���̃`�����N�̖��O���擾���܂��D
    //SizeBuf(ChunkPointer)�́C�`�����N�̐擪�|�C���^����C���̃`�����N��Buffer�̃T�C�Y(size_t)��
    //�擾���܂��D
    //DescendBuf(chunkPointer)�́C�`�����N�̐擪�|�C���^����C���̃`�����N��Buffer�̐擪�|�C���^
    //���擾���܂��D
    //�`�����N���e�ł������ꍇ�́C�e������ԏ��߂̎q��Chunk�̐擪�̃|�C���^�ɂȂ�܂��D
    //NextBuf�́C�`�����N�̐擪�|�C���^����C���̃`�����N�̐擪�|�C���^���擾���܂��D
    
    
    printf("%s\n",NameBuf(dat));
    
    void *tmp = NextBuf(dat);
    printf("%s\n",NameBuf(tmp));
    tmp = DescendBuf(tmp);
    tmp = DescendBuf(tmp);
    printf("%s\n",NameBuf(tmp));
    int nElem = SizeBuf(tmp)/sizeof(int);
    int* tmp2 = (int*)DescendBuf(tmp);
    for (int i=0;i<nElem;i++){
        printf("%d,",tmp2[i]);
    }
    
    tmp = DescendBuf(dat);
    tmp = NextBuf(tmp);
    printf("%s\n",NameBuf(tmp));
    double* tmp3 = (double*)DescendBuf(tmp);
    nElem = SizeBuf(tmp)/sizeof(double);
    for (int i=0;i<nElem;i++){
        printf("%lf,",tmp3[i]);
    }
    
    
    
    return 0;
}