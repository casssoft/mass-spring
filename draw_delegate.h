#ifndef DRAW_DELEGATE_H_
#define DRAW_DELEGATE_H_

#define DDWIDTH 800
#define DDHEIGHT 600


namespace DrawDelegate {

  bool SetupOpenGL();
  void SetLineSize(float size);
  void SetViewMatrix(float*viewmatrix);
  void DrawLines(float* pos, int npos, float* color, int ncolor);
  void BeginFrame();
}
#endif
