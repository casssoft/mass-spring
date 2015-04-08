#ifndef PARTICLE_SYSTEM_H__
#define PARTICLE_SYSTEM_H__

#include "../Eigen/Core"
#include <vector>
class Particle {
 public:
  Eigen::Vector3d x;
  Eigen::Vector3d v;
  Eigen::Vector3d f;
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

class LineSegment {
 public:
 Eigen::Vector3d x1;
 Eigen::Vector3d x2;
};

class ParticleSystem {
 public:
  ParticleSystem();
  void Update(double timestep);
  float* GetPositions2d(int* size);
  void SetupTriangle();
  void SetupBall(double x, double y);
  void SetupMouseSpring(int to);
  void SetMouseSpring(bool enabled);
  void SetMousePos(double x, double y);

  std::vector<Spring> springs;
  std::vector<Particle> particles;
 private:
  void ComputeForces();
  void ExplicitEuler(double timestep);
  std::vector<float> posTemp;
  std::vector<double> phaseTemp;
  int mouseP;
  std::vector<int> mouseSprings;
};
#endif // PARTICLE_SYSTEM_H__

