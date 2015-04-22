#ifndef SCENE_H__
#define SCENE_H__

#include <vector>

class ParticleSystem;
class Scene {
 public:
  Scene();
  void InitTime();
  void EndOfFrame();
  void RecalculateFps();
  double GetTimestep();
  void SetLimitFps(bool enable);
  void ToggleLimitFps();
  double GetFps();
  void DrawScene(ParticleSystem* m, int strainSize, float xpos, float ypos, float zpos);
  void DrawGrid(int gridSize);
 private:
  bool limitFps;
  int GridFloor(double x, int mpG);
  std::vector<float> gridpoints;
  std::vector<float> gridcolors;
  double timestart;
  double curTime;
  int frames;
  double timeAhead;
  double prevTime;
  double secondsPerFrame;
};

#endif
