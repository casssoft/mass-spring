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
  float* GetPositions2d(int* size);
  void SetupTriangle();
  void SetupTriforce();
  void SetupBall(double x, double y);
  void SetupMouseSpring(int to);
  void SetMouseSpring(bool enabled);
  void SetMousePos(double x, double y);
  void Reset();
  void SetSpringProperties(double k, double c);

  std::vector<Spring> springs;
  std::vector<Particle> particles;
 private:
  void ComputeForces();
  void ExplicitEuler(double timestep);
  void ImplicitEuler(double timestep);
  std::vector<float> posTemp;
  std::vector<double> phaseTemp;
  int mouseP;
  std::vector<int> mouseSprings;
  double stiffness;
  double dampness;
};
#endif // PARTICLE_SYSTEM_H__

