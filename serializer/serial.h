#define NAME_STR_SIZE 16

typedef struct{
    char name[NAME_STR_SIZE];
    size_t bufSize;
    int isParent;
} ChunkHeader;

typedef void* pChunkBody;

typedef struct chunkIndex{
    char name[NAME_STR_SIZE];
    size_t bufSize;
    pChunkBody p_body;
    chunkIndex * p_next;
    chunkIndex * p_child;
} ChunkIndex;

ChunkIndex constructChunk(char *name,size_t size,void *p_buf){
    ChunkIndex index;
    strcpy(index.name,name);
    index.bufSize = size;
    index.p_body = p_buf;
    index.p_next = NULL;
    index.p_child = NULL;
    return index;
}

ChunkIndex groupChunk(char *name,ChunkIndex**array,int count){
    ChunkIndex index;
    strcpy(index.name,name);
    index.bufSize = 0;
    index.p_body = NULL;
    index.p_child = array[0];
    index.p_next = NULL;
    
    for (int i =0;i<count-1;i++){
        array[i]->p_next = array[i+1];
    }
    
    return index;
}

int updateChunkInfo(ChunkIndex &index){
    size_t size = 0;
    size += sizeof(ChunkHeader);
    if (index.p_child != NULL){
        size += updateChunkInfo(*index.p_child);
        index.bufSize = size - sizeof(ChunkHeader);
    }
    if (index.p_body != NULL){
        size += index.bufSize;
    }
    if (index.p_next != NULL){
        size += updateChunkInfo(*index.p_next);
    }
    return size;
}

void saveChunk(FILE* fp,ChunkIndex &index){
    ChunkHeader head;
    strcpy(head.name,index.name);
    head.bufSize = index.bufSize;
    head.isParent = (index.p_child != NULL) ? 1 : 0; //1‚È‚çŽq‚ ‚èD
    fwrite(&head,sizeof(ChunkHeader),1,fp);
    
    
    if (index.p_child != NULL){
        saveChunk(fp,*index.p_child);
    }
    if (index.p_body != NULL){
        fwrite(index.p_body,1,index.bufSize,fp);
    }
    if (index.p_next != NULL){
        saveChunk(fp,*index.p_next);
    }
      
}

char* NameBuf(void *buf){
    ChunkHeader * head = (ChunkHeader*)buf;
    return (head->name);
}
size_t SizeBuf(void *buf){
    ChunkHeader * head = (ChunkHeader*)buf;
    return head->bufSize;
}
void* DescendBuf(void* buf){
    return (void*)((char*)buf+sizeof(ChunkHeader));
}
void* NextBuf(void*buf){
    return (void*)((char*)buf+sizeof(ChunkHeader)+SizeBuf(buf));
}