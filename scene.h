#ifndef SCENE_H__
#define SCENE_H__

#include "../Eigen/Core"
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

  void GetCameraRay(double x, double y, Eigen::Vector3d* origin, Eigen::Vector3d* ray);
  bool walkForward, walkBack, walkRight, walkLeft, walkUp, walkDown;
  int drawMode;
  bool slowMode;
  int groundMode;
  bool prevMode;
  bool limitFps;
  std::vector<float> fpsVec;
 private:
  int frames;
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
