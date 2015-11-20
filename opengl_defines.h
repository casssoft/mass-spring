#ifndef OPENGL_DEFINES_H__
#define OPENGL_DEFINES_H__

#ifndef MACOSX
#define glCreateShader DDglCreateShader
#define glShaderSource DDglShaderSource
#define glCompileShader DDglCompileShader
#define glGetShaderiv DDglGetShaderiv
#define glGetShaderInfoLog DDglGetShaderInfoLog
#define glDeleteShader DDglDeleteShader
#define glCreateProgram DDglCreateProgram
#define glAttachShader DDglAttachShader
#define glDetachShader DDglDetachShader
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
#define glDeleteBuffers DDglDeleteBuffers
#define glBufferSubData DDglBufferSubData
#define glUniform2f DDglUniform2f
#define glVertexAttribPointer DDglVertexAttribPointer
#define glUniformMatrix4fv DDglUniformMatrix4fv
#define glUniform3f DDglUniform3f
#define glBlendEquationSeparate DDglBlendEquationSeparate

#ifdef IN_DRAW_DELEGATE_CPP__
PFNGLCREATESHADERPROC glCreateShader = NULL;
PFNGLSHADERSOURCEPROC glShaderSource = NULL;
PFNGLCOMPILESHADERPROC glCompileShader = NULL;
PFNGLGETSHADERIVPROC glGetShaderiv = NULL;
PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog = NULL;
PFNGLDELETESHADERPROC glDeleteShader = NULL;
PFNGLCREATEPROGRAMPROC glCreateProgram = NULL;
PFNGLATTACHSHADERPROC glAttachShader = NULL;
PFNGLDETACHSHADERPROC glDetachShader = NULL;
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
PFNGLDELETEBUFFERSPROC glDeleteBuffers = NULL;
PFNGLBUFFERSUBDATAPROC glBufferSubData = NULL;
PFNGLUNIFORM2FPROC glUniform2f = NULL;
PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer = NULL;
PFNGLUNIFORMMATRIX4FVPROC glUniformMatrix4fv = NULL;
PFNGLUNIFORM3FPROC glUniform3f = NULL;
PFNGLBLENDEQUATIONSEPARATEPROC glBlendEquationSeparate = NULL;
#else
extern PFNGLCREATESHADERPROC glCreateShader;
extern PFNGLSHADERSOURCEPROC glShaderSource;
extern PFNGLCOMPILESHADERPROC glCompileShader ;
extern PFNGLGETSHADERIVPROC glGetShaderiv ;
extern PFNGLGETSHADERINFOLOGPROC glGetShaderInfoLog ;
extern PFNGLDELETESHADERPROC glDeleteShader ;
extern PFNGLCREATEPROGRAMPROC glCreateProgram ;
extern PFNGLATTACHSHADERPROC glAttachShader ;
extern PFNGLDETACHSHADERPROC glDetachShader ;
extern PFNGLLINKPROGRAMPROC glLinkProgram ;
extern PFNGLGETPROGRAMIVPROC glGetProgramiv ;
extern PFNGLGETPROGRAMINFOLOGPROC glGetProgramInfoLog ;
extern PFNGLUSEPROGRAMPROC glUseProgram ;
extern PFNGLDELETEPROGRAMPROC glDeleteProgram ;
extern PFNGLGETUNIFORMLOCATIONPROC glGetUniformLocation ;
extern PFNGLUNIFORM1IPROC glUniform1i ;
extern PFNGLGETATTRIBLOCATIONPROC glGetAttribLocation ;
extern PFNGLENABLEVERTEXATTRIBARRAYPROC glEnableVertexAttribArray ;
extern PFNGLGENBUFFERSPROC glGenBuffers ;
extern PFNGLBINDBUFFERPROC glBindBuffer ;
extern PFNGLBUFFERDATAPROC glBufferData ;
extern PFNGLDELETEBUFFERSPROC glDeleteBuffers ;
extern PFNGLBUFFERSUBDATAPROC glBufferSubData ;
extern PFNGLUNIFORM2FPROC glUniform2f ;
extern PFNGLVERTEXATTRIBPOINTERPROC glVertexAttribPointer ;
extern PFNGLUNIFORMMATRIX4FVPROC glUniformMatrix4fv ;
extern PFNGLUNIFORM3FPROC glUniform3f ;
extern PFNGLBLENDEQUATIONSEPARATEPROC glBlendEquationSeparate ;

#endif // IN_DRAW_DELEGATE_CPP__

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
#endif
