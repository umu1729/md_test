//#pragma comment(linker, "/SUBSYSTEM:WINDOWS /ENTRY:mainCRTStartup")
#include "GL/freeglut.h"
#include <string.h>
#include "cuHeader.h"
#include <stdio.h>

#define WIDTH (640*2)
#define HEIGHT (360*2)

int Width;
int Height;

//平行移動用
float x = 0.0f;
bool flag = false;
//緑
GLfloat green[] = { 0.0, 1.0, 0.0, 1.0 };
//ライトの位置
GLfloat lightpos[] = { -200.0, -150.0, -500.0, 1.0 };


V3Buf v3Buf;

void display(void)
{

 glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
 glViewport(0, 0, Width,Height);
 glMatrixMode(GL_PROJECTION);
 glLoadIdentity();
 //視野角,アスペクト比(ウィンドウの幅/高さ),描画する範囲(最も近い距離,最も遠い距離)
 gluPerspective(60.0, (double)Width / (double)Height, 1.0, 1000.0);
 glMatrixMode(GL_MODELVIEW);
    glLoadIdentity();
 //視点の設定
 gluLookAt(20.0,-20.0,-20.0, //カメラの座標
      20.0,20.0,20.0, // 注視点の座標
     0.0,1.0,0.0); // 画面の上方向を指すベクトル
 //ライトの設定
 glLightfv(GL_LIGHT0, GL_POSITION, lightpos);
 //マテリアルの設定
 glMaterialfv(GL_FRONT, GL_DIFFUSE, green);
 //平行移動
    for (int i=0;i<v3Buf.nElem;i++){
        glPushMatrix();
        glTranslatef((v3Buf.xBuf)[i],(v3Buf.yBuf)[i],(v3Buf.zBuf)[i]);
        glutSolidSphere(0.3,12,12);
        glPopMatrix();
    }
 glutSwapBuffers();
}
void idle(void)
{
 if(flag){x-=0.5f;}else{x+=0.5f;}
 if(x>50.0f)flag=true;
 if(x<-50.0f)flag=false;
 Sleep(1);
 glutPostRedisplay();
}
void Init(){
 glClearColor(0.3f, 0.3f, 0.3f, 1.0f);
 glEnable(GL_DEPTH_TEST);
 glEnable(GL_LIGHTING);
 glEnable(GL_LIGHT0);

}
void NOTHING(){
}
void Loop(V3Buf buf){
    v3Buf = buf;
    
    glutMainLoopEvent();
    display();
    idle();
}
void LoopNoDebug(V3Buf buf){
}
void Resize( int w, int h )
{
	Width = w;
    Height = h;
}
int main(int argc, char **argv)
{
    int handle;
    //strcmp(argv[1],"debug"))
    if (1==1){
        glutInitWindowPosition(100, 100);
        glutInitWindowSize(WIDTH, HEIGHT);
        glutInit(&argc, argv);
        glutInitDisplayMode(GLUT_RGBA | GLUT_DOUBLE);
        handle = glutCreateWindow("球を表示して平行移動");
        glutDisplayFunc(NOTHING);
        glutIdleFunc(NOTHING);
        glutReshapeFunc( Resize );
        Init();
        cuMain(Loop);
    }else{
        cuMain(LoopNoDebug);
    }

    
    
     return 0;
}