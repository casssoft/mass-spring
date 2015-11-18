#ifndef PARTICLE_SYSTEM_H__
#define PARTICLE_SYSTEM_H__

#include "../Eigen/Core"
#include <vector>
class Particle {
 public:
  Eigen::Vector3d x;
  Eigen::Vector3d v;
  Eigen::Vector3d f;
  Eigen::Vector3d lx;
  double iMass;
  bool mark;
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

class CollisionSystem;
class CollisionSystemPQP;
class ParticleSystem {
 public:
  ParticleSystem();
  ~ParticleSystem();
  void Update(double timestep, bool solveWithguess, bool corotational, int groundMode);
  void onMouseDrag(Eigen::Vector3d origin, Eigen::Vector3d ray, double timestep);
  void onMousePress(Eigen::Vector3d origin, Eigen::Vector3d ray);

  float* GetPositions3d(int* size, bool);
  float* GetSurfaceTriangles3d(int* size);
  float* GetAllTriangles3d(int* size);
  float* GetTetCenter(int*size);

  void GetCameraPosAndSize(double* x, double* y, double*z);
  float* GetTriColors(int* size, double strainSize);
  float* GetStrainSurfaceTriColors(int* size, double strainSize);
  float* GetStrainAllTriColors(int* size, double strainSize);
  float* GetColors(int* size, double strainSize, float x, float y, float z);
  float* GetCenterColors(int*size, double strainSize);

  void SetupSingleSpring();
  void SetupBendingBar();
  void SetupArmadillo();
  void SetupMeshFile(const char*filename);
  void Reset();
  void SetSpringProperties(double k, double volumeConservation, double c, double grav, double gStiffness, double mStiffness, bool useRollback);

  void GetProfileInfo(double& triplet, double& fromtriplet, double& solve, double& equationSetupTime);

  std::vector<Tetrahedra> tets;
  std::vector<Particle> particles;
  std::vector<Particle> fixed_points;
  double groundLevel;
 private:
  void HandleCollisions(double timestep);
  void SetupCollisions(double lowestpoint);
  void MakeFixedPoint(int i, std::vector<int>& edges, std::vector<int>& faces);
  void ComputeForces();
  void ExplicitEuler(double timestep);
  void ImplicitEulerSparse(double timestep);

  void CopyIntoStartPos();
  std::vector<Eigen::Vector3d> startPos;
  std::vector<float> posTemp;
  std::vector<float> colorTemp;
  std::vector<double> phaseTemp;
  std::vector<int> faces;
  std::vector<int> facetotet;
  std::vector<int> outsidePoints;

  std::vector<Eigen::Vector3d> prevPos;
  std::vector<Eigen::Vector3d> prevVel;
  std::vector<Eigen::Vector3d> prevFEXT;
#ifdef COLLISION_SELFCCD
  CollisionSystem* colSys;
#endif
#ifdef COLLISION_PQP
  CollisionSystemPQP* colSys;
#endif
  double stiffness;
  double volConserve;
  double dampness;
  double groundStiffness;
  double mouseStiffness;
  double gravity;
  bool useColSys;
  bool corotational;
  bool colRolBack;
  void AddTet(int x1, int x2, int x3, int x4);
  void GetTetP(int i, Particle*& p1, Particle*& p2, Particle*& p3, Particle*& p4);
  void GetPointP(int i, Particle*& x1);
  void CalculateParticleMass(int i, float springMass);
};
#endif // PARTICLE_SYSTEM_H__

