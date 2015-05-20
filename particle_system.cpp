#include <GLFW/glfw3.h>
#include "particle_system.h"
#include "draw_delegate.h"
#include "meshgen.h"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

ParticleSystem::ParticleSystem() {
  stiffness = 100;
  dampness = 10;
  gravity = 9.8;
  ground = false;
}

void ParticleSystem::Update(double timestep, bool implicit, bool solveWithguess) {
  //if (implicit)  ImplicitEuler(timestep);
  //else
  ImplicitEulerSparse(timestep, solveWithguess);

  if (ground) {
    for (int i = 0; i < particles.size(); i++) {
      if (particles[i].x[1] > DDHEIGHT-20) {
        particles[i].x[1] = DDHEIGHT - 20;
        if (particles[i].v[1] > 0) {
            particles[i].v[1] = -.8 * particles[i].v[1];
          /*if (particles[i].v[0] > 0) {
            particles[i].v[0] -= 5000 * timestep;
          } else {
            particles[i].v[0] += 5000 * timestep;
          }*/
        }
      }
    }
  }
}

float* ParticleSystem::GetPositions3d(int* size) {
  *size = springs.size()*6;
  posTemp.resize(*size);
  for (int i = 0; i < springs.size(); i++) {
    Particle* to,*from;
    GetSpringP(i, to, from);

    posTemp[i*6] = ((float)to->x[0]);
    posTemp[i*6 + 1] = ((float)to->x[1]);
    posTemp[i*6 + 2] = ((float)to->x[2]);
    posTemp[i*6 + 3] = ((float)from->x[0]);
    posTemp[i*6 + 4] = ((float)from->x[1]);
    posTemp[i*6 + 5] = ((float)from->x[2]);
  }
  return posTemp.data();
}

void ParticleSystem::GetCameraPosAndSize(double* x, double*y, double* z) {
  *x = 0;
  *y = 0;
  *z = 0;
  for (int i = 0; i < particles.size(); ++i) {
    *x += particles[i].x[0];
    *y += particles[i].x[1];
    *z += particles[i].x[2];
  }
  for (int i = 0; i < fixed_points.size(); ++i) {
    *x += fixed_points[i].x[0];
    *y += fixed_points[i].x[1];
    *z += fixed_points[i].x[2];
  }
  *x /= particles.size() + fixed_points.size();
  *y /= particles.size() + fixed_points.size();
  *z /= particles.size() + fixed_points.size();
}

static void LerpColors(double strain, float*color3) {
   if (strain < 1) {
      *(color3) = 0;
      *(color3+1) = 0;
      *(color3+2) = strain;
   } else if (strain < 2) {
      *(color3) = strain - 1;
      *(color3+1) = 0;
      *(color3+2) = strain;
   } else {
      *(color3) = 1;
      *(color3+1) = 0;
      *(color3+2) = 1 - (strain -2);
   }
}

float* ParticleSystem::GetColors(int* size, int strainSize) {
  *size = springs.size()*6;
  colorTemp.resize(*size);
  for (int i = 0; i < springs.size(); i++) {
    Particle* to,*from;
    GetSpringP(i, to, from);

    float strain = ((to->x - from->x).norm() - springs[i].L) / springs[i].L;
    if (strain < 0) strain *= -1;

    strain *= strainSize;
    LerpColors(strain, &(colorTemp[i*6]));
    LerpColors(strain, &(colorTemp[i*6+3]));
  }
  return colorTemp.data();
}

void ParticleSystem::SetupSingleSpring() {
  Reset();
  particles.emplace_back();
  particles.emplace_back();
  particles[0].x << 0.0, 0.0, 0.0;
  particles[0].v << -1, 0.0, 0.0;
  particles[0].iMass = 1;
  particles[1].x << 0.0, 0.1, 0.0;
  particles[1].v << 1.0, 0.0, 0.0;
  particles[1].iMass = 1;

  springs.emplace_back();
  springs[0].to = 0;
  springs[0].from = 1;
  springs[0].k = stiffness;
  springs[0].L = .1;
  springs[0].c = dampness;
  gravity = 0;
  ground = false;

}

void ParticleSystem::SetupTriangle() {
  Reset();
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles[0].x << 0.0, 0.0, 0.0;
  particles[0].v << 0.0, 0.0, 0.0;
  particles[0].iMass = 1;
  particles[1].x << 1.0, 0.0, 0.0;
  particles[1].v << 0.0, 0.0, 0.0;
  particles[1].iMass = 1;
  particles[2].x << 1.0, 1.0, 0.0;
  particles[2].v << 0.0, 0.0, 0.0;
  particles[2].iMass = 1;

  AddSpring(0, 1);
  AddSpring(0, 2);
  AddSpring(1, 2);
  gravity = 0;
  ground = true;
}

// Assumes this is the first Setup function called
void ParticleSystem::SetupTriforce() {
  Reset();
  // 6 particles
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles[0].x << 20.0, 10.0, 0.0;
  particles[0].v << 0.0, 0.0, 0.0;
  particles[0].iMass = 1;
  particles[1].x << 17.5, 12.5, 0.0;
  particles[1].v << 0.0, 0.0, 0.0;
  particles[1].iMass = 1;
  particles[2].x << 22.5, 12.5, 0.0;
  particles[2].v << 0.0, 0.0, 0.0;
  particles[2].iMass = 1;

  particles[3].x << 15.0, 15.0, 0.0;
  particles[3].v << 0.0, 0.0, 0.0;
  particles[3].iMass = 1;
  particles[4].x << 20.0, 15.0, 0.0;
  particles[4].v << 0.0, 0.0, 0.0;
  particles[4].iMass = 1;
  particles[5].x << 25.0, 15.0, 0.0;
  particles[5].v << 0.0, 0.0, 0.0;
  particles[5].iMass = 1;

  // 9 springs
  springs.emplace_back();
  springs.emplace_back();
  springs.emplace_back();

  springs.emplace_back();
  springs.emplace_back();
  springs.emplace_back();
  springs.emplace_back();
  springs.emplace_back();
  springs.emplace_back();

  springs[0].to = 0;
  springs[0].from = 1;
  springs[0].k = stiffness;
  springs[0].L = 5;
  springs[0].c = dampness;

  springs[1].to = 0;
  springs[1].from = 2;
  springs[1].k = stiffness;
  springs[1].L = 5;
  springs[1].c = dampness;

  springs[2].to = 1;
  springs[2].from = 2;
  springs[2].k = stiffness;
  springs[2].L = 5;
  springs[2].c = dampness;

  springs[3].to = 1;
  springs[3].from = 3;
  springs[3].k = stiffness;
  springs[3].L = 5;
  springs[3].c = dampness;
  springs[4].to = 1;
  springs[4].from = 4;
  springs[4].k = stiffness;
  springs[4].L = 5;
  springs[4].c = dampness;
  springs[5].to = 2;
  springs[5].from = 4;
  springs[5].k = stiffness;
  springs[5].L = 5;
  springs[5].c = dampness;
  springs[6].to = 2;
  springs[6].from = 5;
  springs[6].k = stiffness;
  springs[6].L = 5;
  springs[6].c = dampness;
  springs[7].to = 3;
  springs[7].from = 4;
  springs[7].k = stiffness;
  springs[7].L = 5;
  springs[7].c = dampness;
  springs[8].to = 4;
  springs[8].from = 5;
  springs[8].k = stiffness;
  springs[8].L = 5;
  springs[8].c = dampness;

  fixed_points.emplace_back();
  fixed_points[0].x << 0, 0, 0.0;
  AddSpring(5, -1);
  gravity = 9.8;
  ground = false;
}

void ParticleSystem::SetupBridge(int bridgeL) {
  Reset();
  fixed_points.emplace_back();
  fixed_points.emplace_back();

  //int bridgeL = 10;
  fixed_points[0].x << -4, 0, 0.0;
  fixed_points[1].x << bridgeL*4, 0, 0.0;
  for (int i = 0; i < bridgeL; ++i) {
    particles.emplace_back();
    particles[i].x << i* 4, 0, 0.0;
    particles[i].v << 0, 0, 0.0;
  }
  int l2start = bridgeL;
  for (int i = bridgeL; i < bridgeL*2 +1; i++) {
    particles.emplace_back();
    particles[i].x << (i - bridgeL)* 4 - 2, -4, 0.0;
    particles[i].v << 0, 0, 0.0;
  }

  AddSpring(0, -1);
  AddSpring(bridgeL - 1, -2);
  AddSpring(l2start, -1);
  AddSpring(l2start + bridgeL, -2);

  for (int i = 0; i < bridgeL - 1; ++i) {
    AddSpring(i, i+1);
  }
  for (int i = l2start; i < bridgeL*2; ++i) {
    AddSpring(i, i+1);
  }
  for (int i = 0; i < bridgeL; ++i) {
    AddSpring(l2start + i, i);
    AddSpring(i, l2start + i + 1);
  }
  gravity = 9.8;
  ground = false;
  for (int i = 0; i < particles.size(); ++i) {
    CalculateParticleMass(i, 1);
  }
  /*
  for (int i = 0; i < bridgeL - 1; i++) {
    AddSpring(i*2 + ((i+1)%2), (i+1)*2 + (i%2));
  }

  for (int i = 0; i < bridgeL - 1; i++) {
    AddSpring(i*2, (i+1)*2);
    AddSpring(i*2 + 1, (i+1)*2 + 1);
  }*/
}

void ParticleSystem::SetupBendingBar() {
  Reset();
  int psize;
  double* points;
  std::vector<int> edges;
  MeshGen::GenerateBar(points, psize, edges);


  for (int i = 0; i < psize; ++i) {
    particles.emplace_back();
    particles[i].x << points[i*3], points[i*3 + 1], points[i*3 + 2];
    particles[i].v << 0, 0, 0;
  }
  for (int i = 0; i < psize; ++i) {
    if (particles[i].x[2] == 0) {
      printf("fixed_point!\n");
      MakeFixedPoint(i, edges);
      psize -= 1;
      i--;
    }
  }

  for (int i = 0; i < (edges.size()/2); ++i) {
    AddSpring(edges[i*2], edges[i*2 + 1]);
  }
  for (int i = 0; i < particles.size(); ++i) {
    CalculateParticleMass(i, 200.0/edges.size());
  }
  gravity = 9.8;
  printf("Psize: %d, esize %d\n",psize, edges.size());
  ground = false;
  delete[] points;
}

void ParticleSystem::MakeFixedPoint(int p, std::vector<int>& edges) {
  fixed_points.emplace_back();
  fixed_points[fixed_points.size() - 1].x = particles[p].x;
  for (int i = 0; i < edges.size(); ++i) {
    if (edges[i] == p) {
      edges[i] = -1*fixed_points.size();
    } else if (edges[i] > p) {
      edges[i] -= 1;
    }
  }
  particles.erase(particles.begin() + p);
}

void ParticleSystem::CalculateParticleMass(int i, float springMass) {
  float mass = 0;
  for (int j = 0; j < springs.size(); ++j) {
    if (springs[j].to == i) mass += springMass/2;
    if (springs[j].from == i) mass += springMass/2;
  }
  if (mass == 0) mass = 1;
  particles[i].iMass = 1/mass;
}


void ParticleSystem::SetSpringProperties(double k, double c) {
  stiffness = k;
  dampness = c;
}

void ParticleSystem::ComputeForces() {
  //Zero all forces
  for (int i = 0; i < particles.size(); i++) {
    particles[i].f << 0.0, gravity/particles[i].iMass, 0.0;
  }
  //Compute spring forces
  int sSize = springs.size();
  for (int i = 0; i < sSize; i++) {
    Particle *to, *from;
    GetSpringP(i, to, from);

    Eigen::Vector3d springdir;
    double length = (to->x - from->x).norm();
    if (length != 0) {
      springdir = (to->x - from->x) / length;
    } else {
      springdir << 1, 0;
      length = Eigen::NumTraits<double>::epsilon();
    }
    Eigen::Vector3d force =  (length - springs[i].L) * springs[i].k * springdir;
    // Damping
    force += springs[i].c * springdir * ((to->v - from->v).dot(springdir));
    to->f -= force;
    from->f += force;
  }
}

void ParticleSystem::ExplicitEuler(double timestep) {
  phaseTemp.resize(particles.size() * 6);
  for (int i = 0; i < particles.size(); i++) {
    phaseTemp[i * 6] = particles[i].x[0];
    phaseTemp[i * 6 + 1] = particles[i].x[1];
    phaseTemp[i * 6 + 2] = particles[i].x[2];
    phaseTemp[i * 6 + 3] = particles[i].v[0];
    phaseTemp[i * 6 + 4] = particles[i].v[1];
    phaseTemp[i * 6 + 5] = particles[i].v[2];
  }
  ComputeForces();
  for (int i = 0; i < particles.size(); i++) {
    particles[i].v += particles[i].f * particles[i].iMass * timestep/2;
    particles[i].x += particles[i].v * timestep/2;
  }
  ComputeForces();
  for (int i = 0; i < particles.size(); i++) {
    particles[i].x[0] = phaseTemp[i * 6];
    particles[i].x[1] = phaseTemp[i * 6 + 1];
    particles[i].x[2] = phaseTemp[i * 6 + 2];
    particles[i].v[0] = phaseTemp[i * 6 + 3];
    particles[i].v[1] = phaseTemp[i * 6 + 4];
    particles[i].v[2] = phaseTemp[i * 6 + 5];
    particles[i].v += particles[i].f * particles[i].iMass * timestep;
    particles[i].x += particles[i].v * timestep;
  }
}
namespace {
  void PushbackMatrix3d(std::vector<Eigen::Triplet<double>>& tlist, Eigen::Matrix3d& temp, int startcol, int startrow, int mul) {
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        tlist.push_back(Eigen::Triplet<double>(startcol + i, startrow + j, mul * temp(i, j)));
      }
    }
  }

bool hasPrev = false;
Eigen::VectorXd vdiffprev;

//Eigen::SparseMatrix<double> iesA;
//Eigen::VectorXd iesb;
//
//Eigen::SparseMatrix<double> iesdfdx;
//Eigen::SparseMatrix<double> iesdfdv;
//
//std::vector<Eigen::Triplet<double>> iesdfdxtriplet;
//std::vector<Eigen::Triplet<double>> iesdfdvtriplet;

double curTime;
double tripletTime = 0;
double fromTripletTime = 0;
double solveTime = 0;
};
void ParticleSystem::GetProfileInfo(double& triplet, double& fromTriplet, double& solve) {
   triplet = tripletTime;
   fromTriplet = fromTripletTime;
   solve = solveTime;
}
void ParticleSystem::ImplicitEulerSparse(double timestep, bool solveWithguess) {
  int vSize = 3 * particles.size();
  static Eigen::SparseMatrix<double> iesA;
  static Eigen::VectorXd iesb;

  static Eigen::SparseMatrix<double> iesdfdx;
  static Eigen::SparseMatrix<double> iesdfdv;

  static std::vector<Eigen::Triplet<double>> iesdfdxtriplet;
  static std::vector<Eigen::Triplet<double>> iesdfdvtriplet;

  iesA.resize(vSize, vSize);
  iesb.resize(vSize);
  iesdfdx.resize(vSize, vSize);
  iesdfdv.resize(vSize, vSize);

  Eigen::Matrix3d temp;
  Eigen::Matrix3d tempdv;

  iesdfdxtriplet.clear();
  iesdfdvtriplet.clear();

  curTime = glfwGetTime();
  double tempTime;
  for (int i = 0; i < springs.size(); i++) {
    Particle *to, *from;
    if (springs[i].to < 0)
      to = &(fixed_points[springs[i].to * -1 - 1]);
    else
      to = &(particles[springs[i].to]);
    if (springs[i].from < 0)
      from = &(fixed_points[springs[i].from * -1 - 1]);
    else
      from = &(particles[springs[i].from]);
    Eigen::Vector3d springdir = from->x - to->x;
    double length = springdir.norm();
    if (length == 0)  {
      //printf("zero %d\n", i);
      continue;
    }
    // Jacobian for Hookean spring force
    //temp = ( (springdir * springdir.transpose())/(springdir.transpose() * springdir) + ( Eigen::MatrixXd::Identity(3,3) - (springdir * springdir.transpose())/(springdir.transpose() * springdir)) * ( 1- springs[i].L/length)) * springs[i].k;
    temp = springs[i].k * ( (1 - springs[i].L/length) * (Eigen::MatrixXd::Identity(3,3) - ((springdir/length) * (springdir/length).transpose()))
           + ((springdir/length) * (springdir/length).transpose()));

    tempdv = springs[i].c * ((springdir/length) * (springdir/length).transpose());

    if (springs[i].to >= 0 && springs[i].from >= 0) {
      PushbackMatrix3d(iesdfdxtriplet, temp, springs[i].to * 3, springs[i].from * 3, 1);
      //dfdx.block<3,3>(springs[i].to * 3, springs[i].from * 3) += temp;
      PushbackMatrix3d(iesdfdxtriplet, temp, springs[i].from * 3, springs[i].to * 3, 1);
      //dfdx.block<3,3>(springs[i].from * 3, springs[i].to * 3) += temp;

      PushbackMatrix3d(iesdfdvtriplet, tempdv, springs[i].to * 3, springs[i].from * 3, 1);
      //dfdv.block<3,3>(springs[i].to * 3, springs[i].from * 3) += tempdv;
      PushbackMatrix3d(iesdfdvtriplet, tempdv, springs[i].from * 3, springs[i].to * 3, 1);
      //dfdv.block<3,3>(springs[i].from * 3, springs[i].to * 3) += tempdv;
    }
    if (springs[i].to >= 0) {
      PushbackMatrix3d(iesdfdxtriplet, temp, springs[i].to * 3, springs[i].to * 3, -1);
      //dfdx.block<3,3>(springs[i].to * 3, springs[i].to * 3) -= temp;
      PushbackMatrix3d(iesdfdvtriplet, tempdv, springs[i].to * 3, springs[i].to * 3, -1);
      //dfdv.block<3,3>(springs[i].to * 3, springs[i].to * 3) -= tempdv;
    }
    if (springs[i].from >= 0) {
      PushbackMatrix3d(iesdfdxtriplet, temp, springs[i].from * 3, springs[i].from * 3, -1);
      //dfdx.block<3,3>(springs[i].from * 3, springs[i].from * 3) -= temp;
      PushbackMatrix3d(iesdfdvtriplet, tempdv, springs[i].from * 3, springs[i].from * 3, -1);
      //dfdv.block<3,3>(springs[i].from * 3, springs[i].from * 3) -= tempdv;
    }
  }
  tempTime = glfwGetTime();
  tripletTime += tempTime - curTime;
  curTime = tempTime;
  iesdfdx.setFromTriplets(iesdfdxtriplet.begin(), iesdfdxtriplet.end());
  iesdfdv.setFromTriplets(iesdfdvtriplet.begin(), iesdfdvtriplet.end());

  tempTime = glfwGetTime();
  fromTripletTime += tempTime - curTime;;
  curTime = tempTime;

  ComputeForces();
  Eigen::VectorXd v_0(vSize);
  Eigen::VectorXd f_0(vSize);


  std::vector<Eigen::Triplet<double>> masstriplet;

  for (int i = 0; i < particles.size(); i++) {
    v_0[i * 3] = particles[i].v[0];
    v_0[i * 3 + 1] = particles[i].v[1];
    v_0[i * 3 + 2] = particles[i].v[2];
    f_0[i * 3] = particles[i].f[0];
    f_0[i * 3 + 1] = particles[i].f[1];
    f_0[i * 3 + 2] = particles[i].f[2];
    masstriplet.push_back(Eigen::Triplet<double>(i*3,i*3,1/particles[i].iMass));
    masstriplet.push_back(Eigen::Triplet<double>(i*3+1,i*3+1,1/particles[i].iMass));
    masstriplet.push_back(Eigen::Triplet<double>(i*3+2,i*3+2,1/particles[i].iMass));
    //A.coeffRef(i*3,i*3) = 1/particles[i].iMass;
    //A.coeffRef(i*3+1,i*3+1) = 1/particles[i].iMass;
    //A.coeffRef(i*3+2,i*3+2) = 1/particles[i].iMass;
  }
  iesA.setFromTriplets(masstriplet.begin(), masstriplet.end());
  iesb = timestep * (f_0 + timestep * (iesdfdx * v_0));
  iesA = iesA - (timestep * iesdfdv + timestep * timestep * iesdfdx);
  Eigen::VectorXd vdiff(vSize);
  //Eigen::BiCGSTAB<Eigen::SparseMatrix<double> > cg;
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
  //cg.setMaxIterations(100);
  cg.setTolerance(.001);
  curTime = glfwGetTime();

  cg.compute(iesA);
  if (hasPrev && solveWithguess) vdiff = cg.solveWithGuess(iesb, vdiffprev);
  else vdiff = cg.solve(iesb);

  solveTime += glfwGetTime() - curTime;

  vdiffprev = vdiff;
  hasPrev = true;
  //printf("cg iteratons %d\n", cg.iterations());
  for (int i = 0; i < particles.size(); i++) {
    particles[i].v[0] += vdiff[i * 3];
    particles[i].v[1] += vdiff[i * 3 + 1];
    particles[i].v[2] += vdiff[i * 3 + 2];

    particles[i].x += timestep * particles[i].v;
  }
}

void ParticleSystem::ImplicitEulerSparseLU(double timestep, bool solveWithguess) {
  int vSize = 3 * particles.size();
  static Eigen::SparseMatrix<double> iesA;
  static Eigen::VectorXd iesb;

  static Eigen::SparseMatrix<double> iesdfdx;
  static Eigen::SparseMatrix<double> iesdfdv;

  static std::vector<Eigen::Triplet<double>> iesdfdxtriplet;
  static std::vector<Eigen::Triplet<double>> iesdfdvtriplet;

  iesA.resize(vSize, vSize);
  iesb.resize(vSize);
  iesdfdx.resize(vSize, vSize);
  iesdfdv.resize(vSize, vSize);

  Eigen::Matrix3d temp;
  Eigen::Matrix3d tempdv;

  iesdfdxtriplet.clear();
  iesdfdvtriplet.clear();

  curTime = glfwGetTime();
  double tempTime;
  for (int i = 0; i < springs.size(); i++) {
    Particle *to, *from;
    if (springs[i].to < 0)
      to = &(fixed_points[springs[i].to * -1 - 1]);
    else
      to = &(particles[springs[i].to]);
    if (springs[i].from < 0)
      from = &(fixed_points[springs[i].from * -1 - 1]);
    else
      from = &(particles[springs[i].from]);
    Eigen::Vector3d springdir = from->x - to->x;
    double length = springdir.norm();
    if (length == 0)  {
      //printf("zero %d\n", i);
      continue;
    }
    // Jacobian for Hookean spring force
    //temp = ( (springdir * springdir.transpose())/(springdir.transpose() * springdir) + ( Eigen::MatrixXd::Identity(3,3) - (springdir * springdir.transpose())/(springdir.transpose() * springdir)) * ( 1- springs[i].L/length)) * springs[i].k;
    temp = springs[i].k * ( (1 - springs[i].L/length) * (Eigen::MatrixXd::Identity(3,3) - ((springdir/length) * (springdir/length).transpose()))
           + ((springdir/length) * (springdir/length).transpose()));

    tempdv = springs[i].c * ((springdir/length) * (springdir/length).transpose());

    if (springs[i].to >= 0 && springs[i].from >= 0) {
      PushbackMatrix3d(iesdfdxtriplet, temp, springs[i].to * 3, springs[i].from * 3, 1);
      //dfdx.block<3,3>(springs[i].to * 3, springs[i].from * 3) += temp;
      PushbackMatrix3d(iesdfdxtriplet, temp, springs[i].from * 3, springs[i].to * 3, 1);
      //dfdx.block<3,3>(springs[i].from * 3, springs[i].to * 3) += temp;

      PushbackMatrix3d(iesdfdvtriplet, tempdv, springs[i].to * 3, springs[i].from * 3, 1);
      //dfdv.block<3,3>(springs[i].to * 3, springs[i].from * 3) += tempdv;
      PushbackMatrix3d(iesdfdvtriplet, tempdv, springs[i].from * 3, springs[i].to * 3, 1);
      //dfdv.block<3,3>(springs[i].from * 3, springs[i].to * 3) += tempdv;
    }
    if (springs[i].to >= 0) {
      PushbackMatrix3d(iesdfdxtriplet, temp, springs[i].to * 3, springs[i].to * 3, -1);
      //dfdx.block<3,3>(springs[i].to * 3, springs[i].to * 3) -= temp;
      PushbackMatrix3d(iesdfdvtriplet, tempdv, springs[i].to * 3, springs[i].to * 3, -1);
      //dfdv.block<3,3>(springs[i].to * 3, springs[i].to * 3) -= tempdv;
    }
    if (springs[i].from >= 0) {
      PushbackMatrix3d(iesdfdxtriplet, temp, springs[i].from * 3, springs[i].from * 3, -1);
      //dfdx.block<3,3>(springs[i].from * 3, springs[i].from * 3) -= temp;
      PushbackMatrix3d(iesdfdvtriplet, tempdv, springs[i].from * 3, springs[i].from * 3, -1);
      //dfdv.block<3,3>(springs[i].from * 3, springs[i].from * 3) -= tempdv;
    }
  }
  tempTime = glfwGetTime();
  tripletTime += tempTime - curTime;
  curTime = tempTime;
  iesdfdx.setFromTriplets(iesdfdxtriplet.begin(), iesdfdxtriplet.end());
  iesdfdv.setFromTriplets(iesdfdvtriplet.begin(), iesdfdvtriplet.end());

  tempTime = glfwGetTime();
  fromTripletTime += tempTime - curTime;;
  curTime = tempTime;

  ComputeForces();
  Eigen::VectorXd v_0(vSize);
  Eigen::VectorXd f_0(vSize);


  std::vector<Eigen::Triplet<double>> masstriplet;

  for (int i = 0; i < particles.size(); i++) {
    v_0[i * 3] = particles[i].v[0];
    v_0[i * 3 + 1] = particles[i].v[1];
    v_0[i * 3 + 2] = particles[i].v[2];
    f_0[i * 3] = particles[i].f[0];
    f_0[i * 3 + 1] = particles[i].f[1];
    f_0[i * 3 + 2] = particles[i].f[2];
    masstriplet.push_back(Eigen::Triplet<double>(i*3,i*3,1/particles[i].iMass));
    masstriplet.push_back(Eigen::Triplet<double>(i*3+1,i*3+1,1/particles[i].iMass));
    masstriplet.push_back(Eigen::Triplet<double>(i*3+2,i*3+2,1/particles[i].iMass));
    //A.coeffRef(i*3,i*3) = 1/particles[i].iMass;
    //A.coeffRef(i*3+1,i*3+1) = 1/particles[i].iMass;
    //A.coeffRef(i*3+2,i*3+2) = 1/particles[i].iMass;
  }
  iesA.setFromTriplets(masstriplet.begin(), masstriplet.end());
  iesb = timestep * (f_0 + timestep * (iesdfdx * v_0));
  iesA = iesA - (timestep * iesdfdv + timestep * timestep * iesdfdx);
  Eigen::VectorXd vdiff(vSize);
  //Eigen::SparseLU<Eigen::SparseMatrix<double>, Eigen::COLAMDOrdering<int>> cg;
  Eigen::SimplicialLDLT<Eigen::SparseMatrix<double> > cg;
  curTime = glfwGetTime();

  cg.compute(iesA);
  vdiff = cg.solve(iesb);

  solveTime += glfwGetTime() - curTime;

  for (int i = 0; i < particles.size(); i++) {
    particles[i].v[0] += vdiff[i * 3];
    particles[i].v[1] += vdiff[i * 3 + 1];
    particles[i].v[2] += vdiff[i * 3 + 2];

    particles[i].x += timestep * particles[i].v;
  }
}

/*
void ParticleSystem::ImplicitEulerSolveForNewV(double timestep) {
  int vSize = 2 * particles.size();
  Eigen::MatrixXd A(vSize, vSize);
  Eigen::VectorXd b(vSize);

  Eigen::MatrixXd dfdx(vSize, vSize);
  Eigen::MatrixXd dfdv(vSize, vSize);
  dfdx.setZero();
  dfdv.setZero();
  Eigen::Matrix2d temp;
  Eigen::Matrix2d tempdv;
  for (int i = 0; i < springs.size(); i++) {
    Particle *to, *from;
    if (springs[i].to < 0)
      to = &(fixed_points[springs[i].to * -1 - 1]);
    else
      to = &(particles[springs[i].to]);
    if (springs[i].from < 0)
      from = &(fixed_points[springs[i].from * -1 - 1]);
    else
      from = &(particles[springs[i].from]);
    Eigen::Vector2d springdir = from->x - to->x;
    double length = springdir.norm();
    if (length == 0)  {
      //printf("zero %d\n", i);
      continue;
    }
    // Jacobian for Hookean spring force
    temp = (springs[i].k / (length*length)) * (((length - springs[i].L)/length) * (springdir.dot(springdir)) * Eigen::MatrixXd::Identity(2,2) + (1 - (length - springs[i].L)/length) * springdir * springdir.transpose());
    //temp = springs[i].k * ( (1 - springs[i].L/length) * (Eigen::MatrixXd::Identity(2,2) - ((springdir/length) * (springdir/length).transpose()))
    //       + ((springdir/length) * (springdir/length).transpose()));

    tempdv = springs[i].c * ((springdir/length) * (springdir/length).transpose());

    if (springs[i].to >= 0 && springs[i].from >= 0) {
      dfdx.block<2,2>(springs[i].to * 2, springs[i].from * 2) += temp;
      dfdx.block<2,2>(springs[i].from * 2, springs[i].to * 2) += temp;

      dfdv.block<2,2>(springs[i].to * 2, springs[i].from * 2) += tempdv;
      dfdv.block<2,2>(springs[i].from * 2, springs[i].to * 2) += tempdv;
    }
    if (springs[i].to >= 0) {
      dfdx.block<2,2>(springs[i].to * 2, springs[i].to * 2) -= temp;
      dfdv.block<2,2>(springs[i].to * 2, springs[i].to * 2) -= tempdv;
    }
    if (springs[i].from >= 0) {
      dfdx.block<2,2>(springs[i].from * 2, springs[i].from * 2) -= temp;
      dfdv.block<2,2>(springs[i].from * 2, springs[i].from * 2) -= tempdv;
    }
  }
  ComputeForces();
  Eigen::VectorXd v_0(vSize);
  Eigen::VectorXd f_0(vSize);
  A.setZero();
  for (int i = 0; i < particles.size(); i++) {
    v_0[i * 2] = particles[i].v[0];
    v_0[i * 2 + 1] = particles[i].v[1];
    f_0[i * 2] = particles[i].f[0];
    f_0[i * 2 + 1] = particles[i].f[1];
    A.coeffRef(i*2,i*2) = 1/particles[i].iMass;
    A.coeffRef(i*2+1,i*2+1) = 1/particles[i].iMass;
  }
  b = A * v_0 + timestep * f_0;
  A = A - timestep * timestep * dfdx - timestep *dfdv;
  Eigen::VectorXd vnew(vSize);
  Eigen::ConjugateGradient<Eigen::MatrixXd > cg;
  cg.compute(A);
  vnew = cg.solve(b);
  for (int i = 0; i < particles.size(); i++) {
    particles[i].v[0] = vnew[i * 2];
    particles[i].v[1] = vnew[i * 2 + 1];

    particles[i].x += timestep * particles[i].v;
  }
}*/

void ParticleSystem::AddSpring(int to, int from) {
  springs.emplace_back();
  int index = springs.size() - 1;
  springs[index].to = to;
  springs[index].from = from;
  springs[index].k = stiffness;
  springs[index].c = dampness;
  Particle*s1, *s2;
  GetSpringP(index, s1, s2);
  springs[index].L = (s1->x - s2->x).norm();
}

void ParticleSystem::GetSpringP(int i, Particle*& to, Particle*& from) {
  if (springs[i].to < 0)
    to = &(fixed_points[springs[i].to * -1 - 1]);
  else
    to = &(particles[springs[i].to]);
  if (springs[i].from < 0)
    from = &(fixed_points[springs[i].from * -1 - 1]);
  else
    from = &(particles[springs[i].from]);
}

void ParticleSystem::Reset() {
  particles.clear();
  fixed_points.clear();
  springs.clear();
  hasPrev = false;
}
