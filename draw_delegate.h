#ifndef DRAW_DELEGATE_H_
#define DRAW_DELEGATE_H_

#define DDWIDTH 1500
#define DDHEIGHT 1000


namespace DrawDelegate {

  bool SetupOpenGL();
  void SetLineSize(float size);
  void SetViewMatrix(float*viewmatrix);
  void DrawLines(float* pos, int npos, float* color, int ncolor);
  void DrawTriangles(float* pos, int npos, float* color, int ncolor);
  void DrawPoints(float* pos, int npos, float* color, int ncolor);
  void BeginFrame();
}
#endif
