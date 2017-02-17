@set LIB_GLUT=freeglut\lib\x64
@set CL_DIR="C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin\amd64"
clglut -c gl1.cpp & nvccc -c gl2.cu & nvccc gl1.obj gl2.obj -Xlinker " /LIBPATH:%LIB_GLUT%  " 