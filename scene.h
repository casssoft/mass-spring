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
  void DrawScene(ParticleSystem* m, double strainSize, bool drawPoints);
  void DrawGrid(int gridSize);
  void Update(double timestep);

  bool walkForward, walkBack, walkRight, walkLeft;
  int drawMode;
  bool slowMode;
 private:
  int frames;
  bool limitFps;
  int GridFloor(double x, int mpG);
  std::vector<float> gridpoints;
  std::vector<float> gridcolors;
  double timestart;
  double curTime;
  double timeAhead;
  double prevTime;
  double secondsPerFrame;
  double xpos, ypos, zpos;
  double xtarg, ytarg, ztarg;
};

#endif
