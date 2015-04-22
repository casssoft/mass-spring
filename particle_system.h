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

class ParticleSystem {
 public:
  ParticleSystem();
  void Update(double timestep, bool implicit);
  float* GetPositions3d(int* size);
  void GetCameraPosAndSize(double* x, double* y, double*z);
  float* GetColors(int* size, int strainSize);
  void SetupSingleSpring();
  void SetupTriangle();
  void SetupTriforce();
  void SetupMouseSpring(int to);
  void SetupBridge(int bridgeL);
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
  double stiffness;
  double dampness;
  double gravity;
  bool ground;
  void AddSpring(int to, int from);
  void GetSpringP(int i, Particle*& to, Particle*& from);
  void CalculateParticleMass(int i, float springMass);
};
#endif // PARTICLE_SYSTEM_H__

