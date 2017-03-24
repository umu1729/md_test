#pragma once
#include <stdio.h>
#include <cuda_runtime.h>
#include <time.h>
#include <string.h>

#include "serial.h"


int main(int argc,char** argv){


    //  データを用意する

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
    
    //ChunkIndexの作成
    //ChunkIndexは，対象とするバッファ（配列とか）のポインタや，次のChunkIndex，子のChunkIndexを保持します．
    //ChunkIndexを作成することで，データ構造を定義し，ロード時に使いやすい出力を構成します．
    // constuctChunk(name,size,pointer)は，pointerを先頭として，長さsizeを持つChunkIndexを，
    // nameという名前で作成します．
    // groupChunk(name,ChunkIndexPointerArray,count)は，グループ化したいChunkIndexの，
    // ポインター配列(ChunkIndexPointerArray)と，その配列の長さcountから，nameという名前で
    // グループ化されたChunkIndexを作成します,これは親子構造を構築します．
    
    ChunkIndex chk1 = constructChunk("int1",d1s*sizeof(int),d1);
    ChunkIndex chk2 = constructChunk("flo1",d2s*sizeof(float),d2);
    ChunkIndex* chksA[2] = {&chk1,&chk2};
    ChunkIndex chk3 = groupChunk("g1",chksA,2);
    ChunkIndex chk4 = constructChunk("dou1",d3s*sizeof(double),d3);
    ChunkIndex* chksB[2] = {&chk3,&chk4};
    ChunkIndex chk5 = groupChunk("g2",chksB,2);
    
    //以上の定義で，g2{g1{int1,flo1},dou1}というようなChunkIndexが作成された．
    //ただし，記法　groupName{child1Name,child2Name,...}となっている．
    // ChunkIndexをグループ化することで，出力ファイル内でデータを階層化（親→子→孫）できます．
    
    // Chunk情報の更新
    //updateChunkInfo(topLevelChunkIndex)関数を用いて，ChunkIndexの情報を最新の状態にします．
    //ChunkIndexを新たに作ったり，グループ化した場合は，一番上のChunkIndex(topLevelChunkIndex)を
    //引数に加えて．情報を更新する必要があります．これは内部的には，ファイル出力時のバッファサイズの
    //計算を行っています．
    
    updateChunkInfo(chk5);
    
    //出力用ファイルハンドルを用意して，ChunkIndexを元に，saveChunk(fp,chuk)関数でデータを書き込みます．
    //ChunkIndexがグループ化されたものであれば，グループ全体が書き込まれます．
    //基本的には，一番上の親を引数に渡します
    //（注意：現在では一番上の親が渡されることを想定しているので，子を入れると思わぬバグがあるかもしれません．
    //fpは用意したファイルハンドル，chunkは書き込むChunkIndexです．
    //saveChunkを１度使用するごとにファイルにデータが書き込まれます．
    //このサンプルではチャンクが２回書き込まれています．
    //この場合ではあまり意味がありませんが，例えばsaveChunkの間にChunkIndexによって監視されている
    //ポインターが変更された場合（MDであれば1step進めるなど）
    //saveChunkによって，異なるデータを保存することができます．
    //MDの場合は，chunkindexに位置情報バッファを入れておき，毎stepごとにsaveChunkをすることで，
    //ファイルに各ステップの情報がすべて出力できます．
    
    //出力ファイルの構造は，１つのChunkIndexごとに，Header(データの情報)とBody(データ本体)
    //が連続して書き込まれる構造になっています(Header+BodyでChunkと呼ぶことにします)．
    //親子がない単一のChunkIndexをsaveChunkで何度か書き込む場合は，
    //Header-Body-Header-Body-Header-Body-...という構造になっています．
    //または，Chunk-Chunk-Chunk-...
    //親子関係のあるChunkIndexを保存した場合は，
    //親のBody部分に，子のChunkが埋め込まれます．
    //parent{child}の場合は，
    //Hp-Hc-Bc-Hp-Hc-Bc-... となります Hp:parent's header,Bc:child's header
    //子が複数の場合は，
    //Hp-Hc1-Bc1-Hc2-Bc2-Hp-Hc1-Bc1-Hc2-Bc2-..
    //となります．Headerに，Bodyのサイズの情報が書き込まれているので，簡単な関数を使って
    //子のHeaderや次のHeaderのバッファー位置を特定できます．
    
    FILE * fp= fopen("dat.txt","wb");
    printf("ch:%d\n",(int)sizeof(ChunkHeader));
    saveChunk(fp,chk5);
    saveChunk(fp,chk5);
    fclose(fp);
    
    //ここからは，出力されたファイルから，データを復元するコードです．
    
    printf("LOAD\n");
    
    //まず最初に，出力ファイルをすべて読み込みます．
    //注：すべて読み込まないで徐々に読み込む方式も開発する予定です．
    
    fp = fopen("dat.txt","rb");
    fseek(fp,0,SEEK_END);
    size_t size = ftell(fp);
    fseek(fp,0,SEEK_SET);
    printf("fsize:%d\n",(int)size);
    void *dat = malloc(size);
    fread(dat,1,size,fp);
    
    //読み込むと，バッファポインタが，一番最初に書き込まれたChunkの先頭を指している状態となります．
    //(つまり，Headerの先頭でもある）
    //定義されている便利関数を用いて，バッファポインタを，移動させることによって参照できるようにします．
    //NameBuf(ChunkPointer)は，チャンクの先頭ポインタから，そのチャンクの名前を取得します．
    //SizeBuf(ChunkPointer)は，チャンクの先頭ポインタから，そのチャンクのBufferのサイズ(size_t)を
    //取得します．
    //DescendBuf(chunkPointer)は，チャンクの先頭ポインタから，そのチャンクのBufferの先頭ポインタ
    //を取得します．
    //チャンクが親であった場合は，親が持つ一番初めの子のChunkの先頭のポインタになります．
    //NextBufは，チャンクの先頭ポインタから，次のチャンクの先頭ポインタを取得します．
    
    
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