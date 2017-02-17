typedef struct {
    float* xBuf;
    float* yBuf;
    float* zBuf;
    int nElem;
}V3Buf;
void cuMain( void (*grpc)(V3Buf buf) );