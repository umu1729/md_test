@rem このコンパイルシステムはVC12とWinSDK7.1とfreeglutに依存しています．
@set INC_WSDK="C:\Program Files (x86)\Microsoft SDKs\Windows\v7.1A\Include"
@set INC_GLUT="freeglut\include"
@set INC_VC12="C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\include"
@set LIB_WSDK="C:\Program Files (x86)\Microsoft SDKs\Windows\v7.1A\LIB\x64"
@set LIB_GLUT="freeglut\lib\x64"
@set LIB_VC12="C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\lib\amd64"
@set CL_DIR="C:\Program Files (x86)\Microsoft Visual Studio 12.0\VC\bin\amd64"
%CL_DIR%\cl /I %INC_WSDK% /I %INC_GLUT% /I %INC_VC12% %* /link /LIBPATH:%LIB_WSDK% /LIBPATH:%LIB_GLUT% /LIBPATH:%LIB_VC12% 