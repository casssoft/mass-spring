#include "particle_system.h"
#include "draw_delegate.h"
#include "duckrace.h"
#include "Eigen/Sparse"
#include "Eigen/IterativeLinearSolvers"

ParticleSystem::ParticleSystem() {
  mouseP = -1;
  stiffness = 100;
  dampness = 10;
}

void ParticleSystem::Update(double timestep, bool implicit) {
  if (implicit) ImplicitEuler(timestep);
  else ExplicitEuler(timestep);

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

float* ParticleSystem::GetPositions2d(int* size) {
  *size = springs.size()*4 + 4;
  posTemp.resize(*size);
  for (int i = 0; i < springs.size(); i++) {
    posTemp[i*4] = (float)particles[springs[i].to].x[0];
    posTemp[i*4 + 1] = (float)particles[springs[i].to].x[1];
    posTemp[i*4 + 2] = (float)particles[springs[i].from].x[0];
    posTemp[i*4 + 3] = (float)particles[springs[i].from].x[1];
  }
  posTemp[*size - 4] = 0;
  posTemp[*size - 3] = DDHEIGHT - 20;
  posTemp[*size - 2] = DDWIDTH;
  posTemp[*size - 1] = DDHEIGHT - 20;
  return posTemp.data();
}

// Assumes this is the first Setup function called
void ParticleSystem::SetupTriangle() {
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles[0].x << 100.0, 100.0, 0.0;
  particles[0].v << 0.0, 0.0, 0.0;
  particles[0].iMass = 1;
  particles[1].x << 200.0, 100.0, 0.0;
  particles[1].v << 0.0, 10.0, 0.0;
  particles[1].iMass = 1;
  particles[2].x << 200.0, 200.0, 0.0;
  particles[2].v << 0.0, 0.0, 0.0;
  particles[2].iMass = 1;

  springs.emplace_back();
  springs.emplace_back();
  springs.emplace_back();
  springs[0].to = 0;
  springs[0].from = 1;
  springs[0].k = stiffness;
  springs[0].L = 50;
  springs[0].c = dampness;

  springs[1].to = 0;
  springs[1].from = 2;
  springs[1].k = stiffness;
  springs[1].L = 50;
  springs[1].c = dampness;

  springs[2].to = 1;
  springs[2].from = 2;
  springs[2].k = stiffness;
  springs[2].L = 50;
  springs[2].c = dampness;
}

// Assumes this is the first Setup function called
void ParticleSystem::SetupTriforce() {
  // 6 particles
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  particles[0].x << 200.0, 100.0, 0.0;
  particles[0].v << 0.0, 0.0, 0.0;
  particles[0].iMass = 1;
  particles[1].x << 175.0, 125.0, 0.0;
  particles[1].v << 0.0, 10.0, 0.0;
  particles[1].iMass = 1;
  particles[2].x << 225.0, 125.0, 0.0;
  particles[2].v << 0.0, 0.0, 0.0;
  particles[2].iMass = 1;

  particles[3].x << 150.0, 150.0, 0.0;
  particles[3].v << 0.0, 0.0, 0.0;
  particles[3].iMass = 1;
  particles[4].x << 200.0, 150.0, 0.0;
  particles[4].v << 0.0, 0.0, 0.0;
  particles[4].iMass = 1;
  particles[5].x << 250.0, 150.0, 0.0;
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
  springs[0].L = 50;
  springs[0].c = dampness;

  springs[1].to = 0;
  springs[1].from = 2;
  springs[1].k = stiffness;
  springs[1].L = 50;
  springs[1].c = dampness;

  springs[2].to = 1;
  springs[2].from = 2;
  springs[2].k = stiffness;
  springs[2].L = 50;
  springs[2].c = dampness;

  springs[3].to = 1;
  springs[3].from = 3;
  springs[3].k = stiffness;
  springs[3].L = 50;
  springs[3].c = dampness;
  springs[4].to = 1;
  springs[4].from = 4;
  springs[4].k = stiffness;
  springs[4].L = 50;
  springs[4].c = dampness;
  springs[5].to = 2;
  springs[5].from = 4;
  springs[5].k = stiffness;
  springs[5].L = 50;
  springs[5].c = dampness;
  springs[6].to = 2;
  springs[6].from = 5;
  springs[6].k = stiffness;
  springs[6].L = 50;
  springs[6].c = dampness;
  springs[7].to = 3;
  springs[7].from = 4;
  springs[7].k = stiffness;
  springs[7].L = 50;
  springs[7].c = dampness;
  springs[8].to = 4;
  springs[8].from = 5;
  springs[8].k = stiffness;
  springs[8].L = 50;
  springs[8].c = dampness;
}

void ParticleSystem::SetupBall(double x, double y) {
  Duckrace::MakeBall(this, x, y);
}

void ParticleSystem::SetupMouseSpring(int to) {
  if (mouseP == -1) {
    particles.emplace_back();
    mouseP = particles.size() - 1;
    particles[mouseP].x << 0, 0, 0;
    particles[mouseP].v << 0, 0, 0;
    particles[mouseP].iMass = 1;
  }
  springs.emplace_back();
  mouseSprings.push_back(springs.size() - 1);
  Spring* tempS = &(springs[springs.size() - 1]);
  tempS->to = to;
  tempS->from = mouseP;
  tempS->L = 50;
  tempS->k = stiffness;
  tempS->c = dampness;
}

void ParticleSystem::SetMouseSpring(bool enabled) {
  for (int i = 0; i < mouseSprings.size(); ++i) {
    if (enabled) {
      springs[mouseSprings[i]].k = stiffness;
      springs[mouseSprings[i]].c = dampness;
    } else {
      springs[mouseSprings[i]].k = 0;
      springs[mouseSprings[i]].c = 0;
    }
  }
}


void ParticleSystem::SetMousePos(double x, double y) {
   if (mouseP != -1) {
     particles[mouseP].x << x, y, 0;
     particles[mouseP].v << 0, 0, 0;
   }
}

void ParticleSystem::Reset() {
  mouseP = -1;
  particles.clear();
  springs.clear();
  mouseSprings.clear();
}

void ParticleSystem::SetSpringProperties(double k, double c) {
  stiffness = k;
  dampness = c;
}

void ParticleSystem::ComputeForces() {
  //Zero all forces
  for (int i = 0; i < particles.size(); i++) {
    particles[i].f << 0.0, 200.0, 0.0;
  }
  //Compute spring forces
  int sSize = springs.size();
  for (int i = 0; i < sSize; i++) {
    int to = springs[i].to;
    int from = springs[i].from;
    double length = (particles[to].x - particles[from].x).norm();
    Eigen::Vector3d springdir = (particles[to].x - particles[from].x) / length;
    Eigen::Vector3d force =  (length - springs[i].L) * springs[i].k * springdir;
    // Damping
    force += springs[i].c * springdir * ((particles[to].v - particles[from].v).dot(springdir));
    particles[to].f -= force;
    particles[from].f += force;
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

void ParticleSystem::ImplicitEuler(double timestep) {
  int vSize = 3 * particles.size();
  Eigen::SparseMatrix<double> A(vSize, vSize);
  Eigen::VectorXd b(vSize);

  Eigen::SparseMatrix<double> dfdx(vSize, vSize);
  Eigen::SparseMatrix<double> dfdv(vSize, vSize);
  dfdx.setZero();
  dfdv.setZero();
  Eigen::Matrix3d temp;
  Eigen::Matrix3d tempdv;
  for (int i = 0; i < springs.size(); i++) {
    Eigen::Vector3d springdir = particles[springs[i].from].x - particles[springs[i].to].x;
    double length = springdir.norm();
    // Jacobian for Hookean spring force
    //temp = ( (springdir * springdir.transpose())/(springdir.transpose() * springdir) + ( Eigen::MatrixXd::Identity(3,3) - (springdir * springdir.transpose())/(springdir.transpose() * springdir)) * ( 1- springs[i].L/length)) * springs[i].k;
    temp = springs[i].k * ( (1 - springs[i].L/length) * (Eigen::MatrixXd::Identity(3,3) - ((springdir/length) * (springdir/length).transpose()))
           + ((springdir/length) * (springdir/length).transpose()));

    tempdv = springs[i].c * ((springdir/length) * (springdir/length).transpose());

    for (int i = 0; i < 3; i++) {
      for (int j = 0; j < 3; j++) {
        dfdx.coeffRef(springs[i].to * 3 + i, springs[i].from * 3 + j) += temp(i, j);
        dfdx.coeffRef(springs[i].from * 3 + i, springs[i].to * 3 + j) += temp(i, j);
        dfdx.coeffRef(springs[i].to * 3 + i, springs[i].to * 3 + j) -= temp(i, j);
        dfdx.coeffRef(springs[i].from * 3 + i, springs[i].from * 3 + j) -= temp(i, j);

        dfdv.coeffRef(springs[i].to * 3 + i, springs[i].from * 3 + j) += tempdv(i, j);
        dfdv.coeffRef(springs[i].from * 3 + i, springs[i].to * 3 + j) += tempdv(i, j);
        dfdv.coeffRef(springs[i].to * 3 + i, springs[i].to * 3 + j) -= tempdv(i, j);
        dfdv.coeffRef(springs[i].from * 3 + i, springs[i].from * 3 + j) -= tempdv(i, j);
      }
    }
  }
  ComputeForces();
  Eigen::VectorXd v_0(vSize);
  Eigen::VectorXd f_0(vSize);
  A.setZero();
  for (int i = 0; i < particles.size(); i++) {
    v_0[i * 3] = particles[i].v[0];
    v_0[i * 3 + 1] = particles[i].v[1];
    v_0[i * 3 + 2] = particles[i].v[2];
    f_0[i * 3] = particles[i].f[0];
    f_0[i * 3 + 1] = particles[i].f[1];
    f_0[i * 3 + 2] = particles[i].f[2];
    A.coeffRef(i*3,i*3) = 1/particles[i].iMass;
    A.coeffRef(i*3+1,i*3+1) = 1/particles[i].iMass;
    A.coeffRef(i*3+2,i*3+2) = 1/particles[i].iMass;
  }
  b = timestep * (f_0 + timestep * (dfdx * v_0));
  A -= timestep * dfdv + timestep * timestep * dfdx;
  Eigen::VectorXd vdiff(vSize);
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
  cg.compute(A);
  vdiff = cg.solve(b);
  for (int i = 0; i < particles.size(); i++) {
    particles[i].v[0] += vdiff[i * 3];
    particles[i].v[1] += vdiff[i * 3 + 1];
    particles[i].v[2] += vdiff[i * 3 + 2];

    particles[i].x += timestep * particles[i].v;
  }
}
