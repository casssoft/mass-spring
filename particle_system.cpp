#include "particle_system.h"
#include "draw_delegate.h"
#include "duckrace.h"

ParticleSystem::ParticleSystem() {
  mouseP = -1;
}

void ParticleSystem::Update(double timestep) {
  ExplicitEuler(timestep);
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

  double tempK = 100;
  double tempC = 10;

  springs[0].to = 0;
  springs[0].from = 1;
  springs[0].k = tempK;
  springs[0].L = 50;
  springs[0].c = tempC;

  springs[1].to = 0;
  springs[1].from = 2;
  springs[1].k = tempK;
  springs[1].L = 50;
  springs[1].c = tempC;

  springs[2].to = 1;
  springs[2].from = 2;
  springs[2].k = tempK;
  springs[2].L = 50;
  springs[2].c = tempC;
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
  tempS->k = 100;
  tempS->c = 10;
}

void ParticleSystem::SetMouseSpring(bool enabled) {
  for (int i = 0; i < mouseSprings.size(); ++i) {
    if (enabled) {
      springs[mouseSprings[i]].k = 100;
      springs[mouseSprings[i]].c = 10;
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
    particles[i].x += particles[i].v * timestep/2;
    particles[i].v += particles[i].f * particles[i].iMass * timestep/2;
  }
  ComputeForces();
  for (int i = 0; i < particles.size(); i++) {
    particles[i].x[0] = phaseTemp[i * 6];
    particles[i].x[1] = phaseTemp[i * 6 + 1];
    particles[i].x[2] = phaseTemp[i * 6 + 2];
    particles[i].v[0] = phaseTemp[i * 6 + 3];
    particles[i].v[1] = phaseTemp[i * 6 + 4];
    particles[i].v[2] = phaseTemp[i * 6 + 5];
    particles[i].x += particles[i].v * timestep;
    particles[i].v += particles[i].f * particles[i].iMass * timestep;
  }
}
