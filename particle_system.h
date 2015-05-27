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

class Tetrahedra {
 public:
 int to[4];
 Eigen::Vector3d oldPos[3];
 Eigen::Matrix3d inversePos;
 double posDet;
 double k;
 double c;
 double strain;
};

class ParticleSystem {
 public:
  ParticleSystem();
  void Update(double timestep, bool implicit, bool solveWithguess);
  float* GetPositions3d(int* size);
  float* GetTriangles3d(int* size);
  void GetCameraPosAndSize(double* x, double* y, double*z);
  float* GetTriColors(int* size, int strainSize);
  float* GetColors(int* size, int strainSize);
  void SetupSingleSpring();
  void SetupBendingBar();
  void Reset();
  void SetSpringProperties(double k, double c);

  void GetProfileInfo(double& triplet, double& fromtriplet, double& solve);

  std::vector<Tetrahedra> tets;
  std::vector<Particle> particles;
  std::vector<Particle> fixed_points;
 private:
  void MakeFixedPoint(int i, std::vector<int>& edges);
  void ComputeForces();
  void ExplicitEuler(double timestep);
  void ImplicitEulerSparse(double timestep);

  void CopyIntoStartPos();
  std::vector<Eigen::Vector3d> startPos;
  std::vector<float> posTemp;
  std::vector<float> colorTemp;
  std::vector<double> phaseTemp;
  double stiffness;
  double dampness;
  double gravity;
  bool ground;
  void AddTet(int x1, int x2, int x3, int x4);
  void GetTetP(int i, Particle*& p1, Particle*& p2, Particle*& p3, Particle*& p4);
  void CalculateParticleMass(int i, float springMass);
};
#endif // PARTICLE_SYSTEM_H__

