#ifndef PARTICLE_SYSTEM_H__
#define PARTICLE_SYSTEM_H__

#include "../Eigen/Core"
#include <vector>
class Particle {
 public:
  Eigen::Vector2d x;
  Eigen::Vector2d v;
  Eigen::Vector2d f;
  double iMass;
};

class Spring {
 public:
  int to;
  int from;
  double k; //stiffness
  double L; //rest length
  double c; //damping
};

class ParticleSystem {
 public:
  ParticleSystem();
  void Update(double timestep, bool implicit);
  float* GetPositions2d(int* size, double x, double y, double zoom);
  void GetCameraPosAndSize(double* x, double* y, double* zoom);
  float* GetColors(int* size, int strainSize);
  void SetupSingleSpring();
  void SetupTriangle();
  void SetupTriforce();
  void SetupBall(double x, double y);
  void SetupMouseSpring(int to);
  void SetMouseSpring(bool enabled);
  void SetMousePos(double x, double y);
  void SetupBridge();
  void SetupBridge2(int bridgeL);
  void Reset();
  void SetSpringProperties(double k, double c);

  std::vector<Spring> springs;
  std::vector<Particle> particles;
  std::vector<Particle> fixed_points;
 private:
  void ComputeForces();
  void ExplicitEuler(double timestep);
  void ImplicitEuler(double timestep);
  void ImplicitEulerSolveForNewV(double timestep);
  std::vector<float> posTemp;
  std::vector<float> colorTemp;
  std::vector<double> phaseTemp;
  int mouseP;
  std::vector<int> mouseSprings;
  double stiffness;
  double dampness;
  double gravity;
  bool ground;
  void AddSpring(int to, int from);
  void GetSpringP(int i, Particle*& to, Particle*& from);
};
#endif // PARTICLE_SYSTEM_H__

