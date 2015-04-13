#include "draw_delegate.h"
#include "GLFW/glfw3.h"
#include <stdio.h>
#include <stdlib.h>

#ifndef MACOSX
#define glCreateShader DDglCreateShader
#define glShaderSource DDglShaderSource
#define glCompileShader DDglCompileShader
#define glGetShaderiv DDglGetShaderiv
#define glGetShaderInfoLog DDglGetShaderInfoLog
#define glDeleteShader DDglDeleteShader
#define glCreateProgram DDglCreateProgram
#define glAttachShader DDglAttachShader
#define glLinkProgram DDglLinkProgram
#define glGetProgramiv DDglGetProgramiv
#define glGetProgramInfoLog DDglGetProgramInfoLog
#define glUseProgram DDglUseProgram
#define glDeleteProgram DDglDeleteProgram
#define glGetUniformLocation DDglGetUniformLocation
#define glUniform1i DDglUniform1i
#define glGetAttribLocation DDglGetAttribLocation
#define glEnableVertexAttribArray DDglEnableVertexAttribArray
#define glGenBuffers DDglGenBuffers
#define glBindBuffer DDglBindBuffer
#define glBufferData DDglBufferData
#define glBufferSubData DDglBufferSubData
#define glUniform2f DDglUniform2f
#define glVertexAttribPointer DDglVertexAttribPointer
#define glUniformMatrix4fv DDglUniformMatrix4fv
#define glUniform3f DDglUniform3f

PFNGLCREATESHADERPROC glCreateShader = NULL;
PFNGLSHADERSOURCEPROC glShaderSource = NULL;
PFNGLCOMPILESHADERPROC glCompileShader = NULL;
PFNGLGETSHADERIVPROC glGetShaderiv = NULL;
PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog = NULL;
PFNGLDELETESHADERPROC glDeleteShader = NULL;
PFNGLCREATEPROGRAMPROC glCreateProgram = NULL;
PFNGLATTACHSHADERPROC glAttachShader = NULL;
PFNGLLINKPROGRAMPROC glLinkProgram = NULL;
PFNGLGETPROGRAMIVPROC glGetProgramiv = NULL;
PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog = NULL;
PFNGLUSEPROGRAMPROC glUseProgram = NULL;
PFNGLDELETEPROGRAMPROC glDeleteProgram = NULL;
PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation = NULL;
PFNGLUNIFORM1IPROC glUniform1i = NULL;
PFNGLGETATTRIBLOCATIONPROC glGetAttribLocation = NULL;
PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray = NULL;
PFNGLGENBUFFERSPROC glGenBuffers = NULL;
PFNGLBINDBUFFERPROC glBindBuffer = NULL;
PFNGLBUFFERDATAPROC glBufferData = NULL;
PFNGLBUFFERSUBDATAPROC glBufferSubData = NULL;
PFNGLUNIFORM2FPROC glUniform2f = NULL;
PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer = NULL;
PFNGLUNIFORMMATRIX4FVPROC glUniformMatrix4fv = NULL;
PFNGLUNIFORM3FPROC glUniform3f = NULL;

#define GL_ARRAY_BUFFER                   0x8892
#define GL_STATIC_DRAW                    0x88E4
#define GL_FRAGMENT_SHADER                0x8B30
#define GL_VERTEX_SHADER                  0x8B31
#define GL_COMPILE_STATUS                 0x8B81
#define GL_LINK_STATUS                    0x8B82
#define GL_INFO_LOG_LENGTH                0x8B84
#define GL_TEXTURE0                       0x84C0
#define GL_BGRA                           0x80E1
#define GL_ELEMENT_ARRAY_BUFFER           0x8893
#endif

#define STARTING_SHEETS 16
char * vertexSource =
"#version 120\n"
"attribute vec2 position;\n"
"attribute vec3 colorQ;\n"
"uniform mat4 viewMatrix;\n"
"varying vec3 fragColor;\n"
"void main() {\n"
"  fragColor = colorQ;\n"
"  gl_Position = viewMatrix * (vec4(0.375,0.375,0,0) + vec4(position,1.0,1.0));\n"
"}\n";
char * fragmentSource =
"#version 120\n"
"varying vec3 fragColor;\n"
"void main(void) {\n"
"gl_FragColor = vec4(clamp(fragColor, vec3(0,0,0), vec3(1,1,1)),1.0);\n"
"}\n";

GLuint shaderProgram;
GLuint positionAttribute;
GLuint colorAttribute;
GLuint viewMatrixUniform;
GLuint pointVBO;
GLuint colorVBO;
int height = DDHEIGHT;
int width = DDWIDTH;
float right = width;
float left = 0;
float top = 0;
float bottom = height;
float close = .01;
float away = 100;
float viewmatrix[] = {
  2/(right - left), 0, 0, -1 * (right+left)/(right-left),
  0, 2/(top-bottom), 0, -1 * (top+bottom)/(top-bottom),
  0, 0, -2/(away-close), (away+close)/(away-close),
  0, 0, 0, 1
};

bool DrawDelegate::SetupOpenGL() {
#ifndef MACOSX
        glCreateShader = (PFNGLCREATESHADERPROC) glfwGetProcAddress( "glCreateShader" );
        glShaderSource = (PFNGLSHADERSOURCEPROC) glfwGetProcAddress( "glShaderSource" );
        glCompileShader = (PFNGLCOMPILESHADERPROC) glfwGetProcAddress( "glCompileShader" );
        glGetShaderiv = (PFNGLGETSHADERIVPROC) glfwGetProcAddress( "glGetShaderiv" );
        glGetShaderInfoLog = (PFNGLGETSHADERINFOLOGPROC) glfwGetProcAddress( "glGetShaderInfoLog" );
        glDeleteShader = (PFNGLDELETESHADERPROC) glfwGetProcAddress( "glDeleteShader" );
        glCreateProgram = (PFNGLCREATEPROGRAMPROC) glfwGetProcAddress( "glCreateProgram" );
        glAttachShader = (PFNGLATTACHSHADERPROC) glfwGetProcAddress( "glAttachShader" );
        glLinkProgram = (PFNGLLINKPROGRAMPROC) glfwGetProcAddress( "glLinkProgram" );
        glGetProgramiv = (PFNGLGETPROGRAMIVPROC) glfwGetProcAddress( "glGetProgramiv" );
        glGetProgramInfoLog = (PFNGLGETPROGRAMINFOLOGPROC) glfwGetProcAddress( "glGetProgramInfoLog" );
        glUseProgram = (PFNGLUSEPROGRAMPROC) glfwGetProcAddress( "glUseProgram" );
        glDeleteProgram = (PFNGLDELETEPROGRAMPROC) glfwGetProcAddress( "glDeleteProgram" );
        glGetUniformLocation = (PFNGLGETUNIFORMLOCATIONPROC) glfwGetProcAddress( "glGetUniformLocation" );
        glUniform1i = (PFNGLUNIFORM1IPROC) glfwGetProcAddress( "glUniform1i" );
		glEnableVertexAttribArray = (PFNGLENABLEVERTEXATTRIBARRAYPROC) glfwGetProcAddress( "glEnableVertexAttribArray" );
		glGetAttribLocation = (PFNGLGETATTRIBLOCATIONPROC) glfwGetProcAddress("glGetAttribLocation");
		glGenBuffers = (PFNGLGENBUFFERSPROC) glfwGetProcAddress("glGenBuffers");
		glBindBuffer = (PFNGLBINDBUFFERPROC) glfwGetProcAddress("glBindBuffer");
		glBufferData = (PFNGLBUFFERDATAPROC) glfwGetProcAddress("glBufferData");
		glBufferSubData = (PFNGLBUFFERSUBDATAPROC) glfwGetProcAddress("glBufferSubData");
		glUniform2f = (PFNGLUNIFORM2FPROC) glfwGetProcAddress("glUniform2f");
		glVertexAttribPointer = (PFNGLVERTEXATTRIBPOINTERPROC) glfwGetProcAddress("glVertexAttribPointer");
		glUniformMatrix4fv = (PFNGLUNIFORMMATRIX4FVPROC) glfwGetProcAddress("glUniformMatrix4fv");
		glUniform3f = (PFNGLUNIFORM3FPROC) glfwGetProcAddress("glUniform3f");

		if (!(glActiveTexture && glCreateShader && glShaderSource && glCompileShader &&
		glGetShaderiv && glGetShaderInfoLog && glDeleteShader && glCreateProgram &&
		glAttachShader && glLinkProgram && glGetProgramiv && glGetProgramInfoLog &&
		glUseProgram && glDeleteProgram && glGetUniformLocation && glUniform1i &&
		glGetAttribLocation && glGenBuffers && glBindBuffer && glBufferData && glEnableVertexAttribArray &&
		glUniform2f && glVertexAttribPointer && glUniformMatrix4fv && glUniform3f))
			return 0;
#endif
  //glEnable(GL_TEXTURE_2D);
  //glEnable(GL_ALPHA_TEST);
	//glAlphaFunc(GL_GREATER,0.1f);
  glDisable(GL_DEPTH_TEST);
  glDisable(GL_CULL_FACE);

  /* Set up shaders */
  GLuint vertexShader, fragmentShader;
  vertexShader = glCreateShader(GL_VERTEX_SHADER);
  fragmentShader = glCreateShader(GL_FRAGMENT_SHADER);
  glShaderSource(vertexShader, 1, (const char**)&vertexSource, NULL);
  glShaderSource(fragmentShader, 1,(const char**)&fragmentSource, NULL);
  char buffer[10000];
  int l;
  printf("%s\n%s\n",glGetString(GL_VERSION),glGetString(GL_SHADING_LANGUAGE_VERSION));
  /* Compile our shader objects */
  glCompileShader(vertexShader);
  glCompileShader(fragmentShader);
  glGetShaderInfoLog(vertexShader,10000,&l,buffer);
  if (buffer[0] != 0)  {
    //MessageBox(NULL,buffer,"VERTEX", MB_OK);
    fprintf(stderr, "VERTEX:\n%s",buffer);
  }
  glGetShaderInfoLog(fragmentShader,10000,&l,buffer);
  if (buffer[0] != 0) {
    //MessageBox(NULL,buffer,"FRAGMENT", MB_OK);
    fprintf(stderr, "FRAGMENT:\n%s",buffer);
  }
  shaderProgram = glCreateProgram();
  /* Attach our shaders to our program */
  glAttachShader(shaderProgram, vertexShader);
  glAttachShader(shaderProgram, fragmentShader);
  /* Bind attribute index 0 (shaderAtribute) to in_Position*/
  /* "in_Position" will represent "data" array's contents in the vertex shader */
  //glBindAttribLocation(shaderProgram, shaderAtribute, "in_Position");
  /* Link shader program*/
  glLinkProgram(shaderProgram);
  GLint status;
  glGetProgramiv (shaderProgram, GL_LINK_STATUS, &status);
  if (status == GL_FALSE)
  {
    GLint infoLogLength;
    glGetProgramiv(shaderProgram, GL_INFO_LOG_LENGTH, &infoLogLength);
    char *strInfoLog = new char[infoLogLength + 1];
    glGetProgramInfoLog(shaderProgram, infoLogLength, NULL, strInfoLog);
    fprintf(stderr, "Linker failure: %s\n", strInfoLog);
    delete[] strInfoLog;

    exit(1);
  }
  /* Set shader program as being actively used */
  glUseProgram(shaderProgram);


  positionAttribute = glGetAttribLocation(shaderProgram,"position");
  colorAttribute = glGetAttribLocation(shaderProgram,"colorQ");
  viewMatrixUniform = glGetUniformLocation(shaderProgram, "viewMatrix");
  glEnableVertexAttribArray(positionAttribute);
  glEnableVertexAttribArray(colorAttribute);
  /*-------------------------------------------------------------------------------------------------------*/
  glClearColor(1,1,1,1);


  glGenBuffers(1, &colorVBO);
  glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*1000, NULL, GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glGenBuffers(1, &pointVBO);
  glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*900, NULL, GL_DYNAMIC_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glClearColor(1,1,1,1);
  glLineWidth(5);
  return 1;
}

void DrawDelegate::SetLineSize(float size) {
  glLineWidth(size);
}

void DrawDelegate::BeginFrame() {
  glClear(GL_COLOR_BUFFER_BIT);
  glColor4f(1.0,1.0,1.0,1.0);
  glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
  glVertexAttribPointer(colorAttribute, 3, GL_FLOAT, false, 0, 0);
  glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
  glVertexAttribPointer(positionAttribute, 2, GL_FLOAT, false, 0, 0);
  glUniformMatrix4fv(viewMatrixUniform, 1, true, viewmatrix);
}

namespace {
  int DDbufferSize = 1000;
  int DDcolorbufferSize = 900;
};
void DrawDelegate::DrawLines(float* pos, int npos, float* color, int ncolor) {
  while(npos > DDbufferSize)
    DDbufferSize *= 2;
  while (ncolor > DDcolorbufferSize)
    DDcolorbufferSize *= 2;
  glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDcolorbufferSize, NULL, GL_DYNAMIC_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * ncolor, color);
  glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDbufferSize, NULL, GL_DYNAMIC_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * npos, pos);

  glDrawArrays(GL_LINES, 0, npos/2);
  GLenum temp = glGetError();
  if (temp != GL_NO_ERROR) {
    fprintf(stderr, "Got gl error: %d\n", temp);
  }
}
