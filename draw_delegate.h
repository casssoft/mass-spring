#ifndef DRAW_DELEGATE_H_
#define DRAW_DELEGATE_H_

#define DDWIDTH 800
#define DDHEIGHT 600


namespace DrawDelegate {

  bool SetupOpenGL();
  void DrawLines(float* pos, int npos);
  void BeginFrame();
}
#endif
