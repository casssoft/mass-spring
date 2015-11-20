#include "draw_delegate.h"
#include "GLFW/glfw3.h"
#include <stdio.h>
#include <stdlib.h>
#define IN_DRAW_DELEGATE_CPP__
#include "opengl_defines.h"

#define STARTING_SHEETS 16
char * vertexSource =
"#version 120\n"
"attribute vec3 position;\n"
"attribute vec3 colorQ;\n"
"uniform mat4 viewMatrix;\n"
"varying vec3 fragColor;\n"
"void main() {\n"
"  fragColor = colorQ;\n"
"  gl_Position = viewMatrix * (vec4(position,1.0));\n"
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

namespace {
  int DDbufferSize = 10000;
  int DDcolorbufferSize = 10000;
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
        glDetachShader = (PFNGLDETACHSHADERPROC) glfwGetProcAddress( "glDetachShader" );
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
		glDeleteBuffers = (PFNGLDELETEBUFFERSPROC) glfwGetProcAddress("glDeleteBuffers");
		glBufferSubData = (PFNGLBUFFERSUBDATAPROC) glfwGetProcAddress("glBufferSubData");
		glUniform2f = (PFNGLUNIFORM2FPROC) glfwGetProcAddress("glUniform2f");
		glVertexAttribPointer = (PFNGLVERTEXATTRIBPOINTERPROC) glfwGetProcAddress("glVertexAttribPointer");
		glUniformMatrix4fv = (PFNGLUNIFORMMATRIX4FVPROC) glfwGetProcAddress("glUniformMatrix4fv");
		glUniform3f = (PFNGLUNIFORM3FPROC) glfwGetProcAddress("glUniform3f");
		glBlendEquationSeparate = (PFNGLBLENDEQUATIONSEPARATEPROC) glfwGetProcAddress("glBlendEquationSeparate");

		if (!(glActiveTexture && glCreateShader && glShaderSource && glCompileShader &&
		glGetShaderiv && glGetShaderInfoLog && glDeleteShader && glCreateProgram &&
		glAttachShader &&  glDetachShader && glLinkProgram && glGetProgramiv && glGetProgramInfoLog &&
		glUseProgram && glDeleteProgram && glGetUniformLocation && glUniform1i &&
		glGetAttribLocation && glGenBuffers && glBindBuffer && glBufferData && glDeleteBuffers && glEnableVertexAttribArray &&
		glUniform2f && glVertexAttribPointer && glUniformMatrix4fv && glUniform3f))
			return 0;
#endif
  //glEnable(GL_TEXTURE_2D);
  //glEnable(GL_ALPHA_TEST);
	//glAlphaFunc(GL_GREATER,0.1f);
  glEnable(GL_DEPTH_TEST);
  //glDisable(GL_CULL_FACE);

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


  glGenBuffers(1, &colorVBO);
  glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDcolorbufferSize, NULL, GL_STREAM_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glGenBuffers(1, &pointVBO);
  glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDbufferSize, NULL, GL_STREAM_DRAW);
  glBindBuffer(GL_ARRAY_BUFFER, 0);

  glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
  glVertexAttribPointer(colorAttribute, 3, GL_FLOAT, false, 0, 0);
  glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
  glVertexAttribPointer(positionAttribute, 3, GL_FLOAT, false, 0, 0);
  glClearColor(1,1,1,1);
  glColor4f(1.0,1.0,1.0,1.0);
  //glLineWidth(5);
  glPointSize(10);
  GLenum temp = glGetError();
  if (temp != GL_NO_ERROR) {
    fprintf(stderr, "Got GL error: %x line %i\n", temp, __LINE__);
  }
  return 1;
}

void DrawDelegate::SetLineSize(float size) {
  glLineWidth(size);
}

void DrawDelegate::SetViewMatrix(float* viewmatrix) {
  glUniformMatrix4fv(viewMatrixUniform, 1, false, viewmatrix);
}
void DrawDelegate::BeginFrame() {
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);
  glUseProgram(shaderProgram);

  glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
  glVertexAttribPointer(colorAttribute, 3, GL_FLOAT, false, 0, 0);
  glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
  glVertexAttribPointer(positionAttribute, 3, GL_FLOAT, false, 0, 0);
  glClearColor(1,1,1,1);
  glColor4f(1.0,1.0,1.0,1.0);
}

void DrawDelegate::DrawTriangles(float* pos, int npos, float* color, int ncolor) {
  while(npos > DDbufferSize) {
   printf("buffer not big enough\n");
   DDbufferSize *= 2;
  }
  while (ncolor > DDcolorbufferSize)
    DDcolorbufferSize *= 2;
  glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDcolorbufferSize, NULL, GL_STREAM_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * ncolor, color);
  glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDbufferSize, NULL, GL_STREAM_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * npos, pos);

  glDrawArrays(GL_TRIANGLES, 0, npos/3);
  GLenum temp = glGetError();
  if (temp != GL_NO_ERROR) {
    fprintf(stderr, "Got GL error: %x line %i\n", temp, __LINE__);
  }
}

void DrawDelegate::DrawLines(float* pos, int npos, float* color, int ncolor) {
  while(npos > DDbufferSize) {
   printf("buffer not big enough\n");
   DDbufferSize *= 2;
  }
  while (ncolor > DDcolorbufferSize)
    DDcolorbufferSize *= 2;
  glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDcolorbufferSize, NULL, GL_STREAM_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * ncolor, color);
  glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDbufferSize, NULL, GL_STREAM_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * npos, pos);

  glDrawArrays(GL_LINES, 0, npos/3);
  GLenum temp = glGetError();
  if (temp != GL_NO_ERROR) {
    fprintf(stderr, "Got GL error: %x line %i\n", temp, __LINE__);
  }
}

void DrawDelegate::DrawPoints(float* pos, int npos, float* color, int ncolor) {
  while(npos > DDbufferSize) {
   printf("buffer not big enough\n");
   DDbufferSize *= 2;
  }
  while (ncolor > DDcolorbufferSize)
    DDcolorbufferSize *= 2;
  glBindBuffer(GL_ARRAY_BUFFER, colorVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDcolorbufferSize, NULL, GL_STREAM_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * ncolor, color);
  glBindBuffer(GL_ARRAY_BUFFER, pointVBO);
  glBufferData(GL_ARRAY_BUFFER, sizeof(float)*DDbufferSize, NULL, GL_STREAM_DRAW);
  glBufferSubData(GL_ARRAY_BUFFER, 0, sizeof(float) * npos, pos);

  glDrawArrays(GL_POINTS, 0, npos/3);
  GLenum temp = glGetError();
  if (temp != GL_NO_ERROR) {
    //fprintf(stderr, "%s\n", gluErrorString(error));
    fprintf(stderr, "Got GL error: %x line %i\n", temp, __LINE__);
  }
}
