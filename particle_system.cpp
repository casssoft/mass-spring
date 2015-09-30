#include <GLFW/glfw3.h>
#include "particle_system.h"
#include "draw_delegate.h"
#include "meshgen.h"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include <iostream>

ParticleSystem::ParticleSystem() {
  stiffness = 1000;
  dampness = 10;
  gravity = 9.8;
  ground = false;
}

void ParticleSystem::Update(double timestep, bool solveWithguess, bool coro) {
  // Always solveWithguess
  // coro means whether to use corotational linear FEM or normal linear FEM
  corotational = coro;
  ImplicitEulerSparse(timestep);
  //ExplicitEuler(timestep);

  // Optionally make things bounce of the ground
  if (ground) {
    for (int i = 0; i < particles.size(); i++) {
      if (particles[i].x[1] > 5) {
        particles[i].x[1] = 5;
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

// Help function to show strain properly through color
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

// Get the center of all tets
float* ParticleSystem::GetTetCenter(int* size) {
  *size = tets.size() * 3;
  posTemp.resize(*size);
  for (int i = 0; i < *size/3; i++) {
    Particle *p1, *p2, *p3, *p4;
    GetTetP(i, p1, p2, p3, p4);
    Eigen::Vector3d avg = (p1->x + p2->x + p3->x + p4->x)/4;
    posTemp[i*3] = avg[0];
    posTemp[i*3+1] = avg[1];
    posTemp[i*3+2] = avg[2];
  }
  return posTemp.data();
}

// Get the strain colors of all tets
float* ParticleSystem::GetCenterColors(int*size, double strainSize) {
  *size = tets.size() * 3;
  colorTemp.resize(*size);
  for (int i = 0; i < *size/3; i++) {
    Particle *p1, *p2, *p3, *p4;
    GetTetP(i, p1, p2, p3, p4);
    Eigen::Matrix3d temp;
    temp << p2->x - p1->x, p3->x - p1->x, p4->x - p1->x;
    Eigen::Matrix3d deformGradient = (temp * tets[i].inversePos) - Eigen::Matrix3d::Identity();
    Eigen::Matrix3d greenStrainTensor;
    if (corotational) {
      greenStrainTensor = .5 * (deformGradient + deformGradient.transpose() +
                                                deformGradient.transpose() * deformGradient);
    } else {
      greenStrainTensor = .5 * (deformGradient + deformGradient.transpose());
    }
    double v = .4;
    Eigen::VectorXd strainVec(6);
    strainVec << greenStrainTensor(0,0), greenStrainTensor(1,1),
                                    greenStrainTensor(2,2), greenStrainTensor(1,0),
                                    greenStrainTensor(1,2), greenStrainTensor(2,0);

    Eigen::MatrixXd strainToStress(6,6);
    strainToStress << 1 - v, v, v, 0, 0, 0,
                                           v, 1 - v, v, 0, 0, 0,
                                           v, v, 1 - v, 0, 0, 0,
                                           0, 0, 0, 1 - 2*v, 0, 0,
                                           0, 0, 0, 0, 1 - 2*v, 0,
                                           0, 0, 0, 0, 0, 1 - 2*v;
    Eigen::VectorXd stressVec(6);
    stressVec = (tets[i].k/((1 + v) * (1 - 2*v))) * strainToStress * strainVec;
    Eigen::Matrix3d stressTensor;
    stressTensor << stressVec[0], stressVec[3], stressVec[5],
                    stressVec[3], stressVec[1], stressVec[4],
                    stressVec[5], stressVec[4], stressVec[2];
    double strain = stressTensor.norm();

    LerpColors(strain * strainSize, &(colorTemp[i*3]));
  }
  return colorTemp.data();
}

// Get colors for lines
float* ParticleSystem::GetColors(int* size, double strainSize, float xpos, float ypos, float zpos) {
  *size = tets.size()* 6 * 3 * 2;
  int perTet = 6 * 3 * 2;
  colorTemp.resize(*size);
  Eigen::Vector3d campos;
  campos << xpos, ypos, zpos;
  for (int i = 0; i < *size/3; i++) {
    Eigen::Vector3d ppos;
    ppos << posTemp[i*3], posTemp[i*3+1], posTemp[i*3+2];
    double color = 1/(ppos - campos).norm();

    colorTemp[i*3] = color + (i%2)/2.0;
    colorTemp[i*3 + 1] = color;
    colorTemp[i*3 + 2] = color + ((i+1)%2)/2.0;

    //colorTemp[i*6 + 3] = 1;
    //colorTemp[i*6 + 4] = 0;
    //colorTemp[i*6 + 5] = 0;
  }
  return colorTemp.data();
}

// Get positions for lines
float* ParticleSystem::GetPositions3d(int* size) {
  *size = tets.size()* 6 * 3 * 2;
  int perTet = 6 * 3 * 2;
  posTemp.resize(*size);
  for (int i = 0; i < tets.size(); i++) {
    Particle *p1, *p2, *p3, *p4;
    GetTetP(i, p1, p2, p3, p4);
    int c = 0;
    // p1 to p2
    posTemp[i*perTet + c++] = ((float)p1->x[0]);
    posTemp[i*perTet + c++] = ((float)p1->x[1]);
    posTemp[i*perTet + c++] = ((float)p1->x[2]);

    posTemp[i*perTet + c++] = ((float)p2->x[0]);
    posTemp[i*perTet + c++] = ((float)p2->x[1]);
    posTemp[i*perTet + c++] = ((float)p2->x[2]);

    // p1 to p3
    posTemp[i*perTet + c++] = ((float)p1->x[0]);
    posTemp[i*perTet + c++] = ((float)p1->x[1]);
    posTemp[i*perTet + c++] = ((float)p1->x[2]);

    posTemp[i*perTet + c++] = ((float)p3->x[0]);
    posTemp[i*perTet + c++] = ((float)p3->x[1]);
    posTemp[i*perTet + c++] = ((float)p3->x[2]);

    // p1 to p4
    posTemp[i*perTet + c++] = ((float)p1->x[0]);
    posTemp[i*perTet + c++] = ((float)p1->x[1]);
    posTemp[i*perTet + c++] = ((float)p1->x[2]);

    posTemp[i*perTet + c++] = ((float)p4->x[0]);
    posTemp[i*perTet + c++] = ((float)p4->x[1]);
    posTemp[i*perTet + c++] = ((float)p4->x[2]);

    // p2 to p3
    posTemp[i*perTet + c++] = ((float)p2->x[0]);
    posTemp[i*perTet + c++] = ((float)p2->x[1]);
    posTemp[i*perTet + c++] = ((float)p2->x[2]);

    posTemp[i*perTet + c++] = ((float)p3->x[0]);
    posTemp[i*perTet + c++] = ((float)p3->x[1]);
    posTemp[i*perTet + c++] = ((float)p3->x[2]);

    // p2 to p4
    posTemp[i*perTet + c++] = ((float)p2->x[0]);
    posTemp[i*perTet + c++] = ((float)p2->x[1]);
    posTemp[i*perTet + c++] = ((float)p2->x[2]);

    posTemp[i*perTet + c++] = ((float)p4->x[0]);
    posTemp[i*perTet + c++] = ((float)p4->x[1]);
    posTemp[i*perTet + c++] = ((float)p4->x[2]);

    // p3 to p4
    posTemp[i*perTet + c++] = ((float)p3->x[0]);
    posTemp[i*perTet + c++] = ((float)p3->x[1]);
    posTemp[i*perTet + c++] = ((float)p3->x[2]);

    posTemp[i*perTet + c++] = ((float)p4->x[0]);
    posTemp[i*perTet + c++] = ((float)p4->x[1]);
    posTemp[i*perTet + c++] = ((float)p4->x[2]);
  }
  return posTemp.data();
}

// Get triangler colors, its just a light calculation
float* ParticleSystem::GetTriColors(int* size, double strainSize) {
  *size = faces.size()*3;
  colorTemp.resize(*size);
  int perTri = 9;
  Eigen::Vector3d light1,light2;
  light1 << 1, 1, -1;
  light1.normalize();
  light2 << -1, 1, -1;
  light2.normalize();
  for (int i = 0; i < *size/9; i++) {
    int c = 0;
    Eigen::Vector3d p0, p1, p2;
    p0 << posTemp[i*perTri + c++], posTemp[i*perTri + c++], posTemp[i*perTri + c++];
    p1 << posTemp[i*perTri + c++], posTemp[i*perTri + c++], posTemp[i*perTri + c++];
    p2 << posTemp[i*perTri + c++], posTemp[i*perTri + c++], posTemp[i*perTri + c++];
    p1 -= p0;
    p2 -= p0;
    p0 = p1.cross(p2);
    p0.normalize();
    double color1 = p0.dot(light1);
    double color2 = 0;//.5*p0.dot(light2);
    int c2 = 0;
    colorTemp[i*perTri + c2++] = color1+color2;
    colorTemp[i*perTri + c2++] = color1+color2;
    colorTemp[i*perTri + c2++] = color1+color2;

    colorTemp[i*perTri + c2++] = color1+color2;
    colorTemp[i*perTri + c2++] = color1+color2;
    colorTemp[i*perTri + c2++] = color1+color2;

    colorTemp[i*perTri + c2++] = color1+color2;
    colorTemp[i*perTri + c2++] = color1+color2;
    colorTemp[i*perTri + c2++] = color1+color2;
    //LerpColors((i/(4*3))/10, &(colorTemp[i*3])); 
  }
  return colorTemp.data();
}

// Get the strain for the surfarce of tris
float* ParticleSystem::GetStrainSurfaceTriColors(int* size, double strainSize) {
  *size = faces.size()*3;
  colorTemp.resize(*size);
  for (int i = 0; i < tets.size(); i++) {
    Particle *p1,*p2,*p3,*p4;
    GetTetP(i, p1, p2, p3, p4);

    Eigen::Matrix3d temp;
    temp << p2->x - p1->x, p3->x - p1->x, p4->x - p1->x;
    Eigen::Matrix3d deformGradient = (temp * tets[i].inversePos) - Eigen::Matrix3d::Identity();
    Eigen::Matrix3d greenStrainTensor;
    if (corotational) {
      greenStrainTensor = .5 * (deformGradient + deformGradient.transpose() +
                                                deformGradient.transpose() * deformGradient);
    } else {
      greenStrainTensor = .5 * (deformGradient + deformGradient.transpose());
    }
    double v = .4;
    Eigen::VectorXd strainVec(6);
    strainVec << greenStrainTensor(0,0), greenStrainTensor(1,1),
                                    greenStrainTensor(2,2), greenStrainTensor(1,0),
                                    greenStrainTensor(1,2), greenStrainTensor(2,0);

    Eigen::MatrixXd strainToStress(6,6);
    strainToStress << 1 - v, v, v, 0, 0, 0,
                                           v, 1 - v, v, 0, 0, 0,
                                           v, v, 1 - v, 0, 0, 0,
                                           0, 0, 0, 1 - 2*v, 0, 0,
                                           0, 0, 0, 0, 1 - 2*v, 0,
                                           0, 0, 0, 0, 0, 1 - 2*v;
    Eigen::VectorXd stressVec(6);
    stressVec = (tets[i].k/((1 + v) * (1 - 2*v))) * strainToStress * strainVec;
    Eigen::Matrix3d stressTensor;
    stressTensor << stressVec[0], stressVec[3], stressVec[5],
                    stressVec[3], stressVec[1], stressVec[4],
                    stressVec[5], stressVec[4], stressVec[2];
    double strain = stressTensor.norm();
    tets[i].strain = strain;
  }
  for (int i = 0; i < faces.size()/3; i++) {
    double strain = tets[facetotet[i]].strain;
    LerpColors(strain * strainSize, &(colorTemp[i*9]));
    LerpColors(strain * strainSize, &(colorTemp[i*9+3]));
    LerpColors(strain * strainSize, &(colorTemp[i*9+6]));
  }
  return colorTemp.data();
}

float* ParticleSystem::GetSurfaceTriangles3d(int* size) {
  *size = faces.size()*3;
  posTemp.resize(*size);
  for (int i = 0; i < faces.size(); i++) {
    Particle* p1;
    GetPointP(faces[i], p1);
    posTemp[i*3] = p1->x[0];
    posTemp[i*3+1] = p1->x[1];
    posTemp[i*3+2] = p1->x[2];
 }
 return posTemp.data();
}

float* ParticleSystem::GetAllTriangles3d(int* size) {
  *size = tets.size()* 4 * 3 * 3;
  int perTet = 4 * 3 * 3;
  posTemp.resize(*size);
  for (int i = 0; i < tets.size(); i++) {
    Particle *p1, *p2, *p3, *p4;
    GetTetP(i, p1, p2, p3, p4);
    int c = 0;
    // p1 p2 p3
    posTemp[i*perTet + c++] = ((float)p1->x[0]);
    posTemp[i*perTet + c++] = ((float)p1->x[1]);
    posTemp[i*perTet + c++] = ((float)p1->x[2]);

    posTemp[i*perTet + c++] = ((float)p2->x[0]);
    posTemp[i*perTet + c++] = ((float)p2->x[1]);
    posTemp[i*perTet + c++] = ((float)p2->x[2]);

    posTemp[i*perTet + c++] = ((float)p3->x[0]);
    posTemp[i*perTet + c++] = ((float)p3->x[1]);
    posTemp[i*perTet + c++] = ((float)p3->x[2]);

    // p1 p4 p2
    posTemp[i*perTet + c++] = ((float)p1->x[0]);
    posTemp[i*perTet + c++] = ((float)p1->x[1]);
    posTemp[i*perTet + c++] = ((float)p1->x[2]);

    posTemp[i*perTet + c++] = ((float)p4->x[0]);
    posTemp[i*perTet + c++] = ((float)p4->x[1]);
    posTemp[i*perTet + c++] = ((float)p4->x[2]);

    posTemp[i*perTet + c++] = ((float)p2->x[0]);
    posTemp[i*perTet + c++] = ((float)p2->x[1]);
    posTemp[i*perTet + c++] = ((float)p2->x[2]);

    // p1 p3 p4
    posTemp[i*perTet + c++] = ((float)p1->x[0]);
    posTemp[i*perTet + c++] = ((float)p1->x[1]);
    posTemp[i*perTet + c++] = ((float)p1->x[2]);

    posTemp[i*perTet + c++] = ((float)p3->x[0]);
    posTemp[i*perTet + c++] = ((float)p3->x[1]);
    posTemp[i*perTet + c++] = ((float)p3->x[2]);

    posTemp[i*perTet + c++] = ((float)p4->x[0]);
    posTemp[i*perTet + c++] = ((float)p4->x[1]);
    posTemp[i*perTet + c++] = ((float)p4->x[2]);

    // p3 p2 p4
    posTemp[i*perTet + c++] = ((float)p3->x[0]);
    posTemp[i*perTet + c++] = ((float)p3->x[1]);
    posTemp[i*perTet + c++] = ((float)p3->x[2]);

    posTemp[i*perTet + c++] = ((float)p2->x[0]);
    posTemp[i*perTet + c++] = ((float)p2->x[1]);
    posTemp[i*perTet + c++] = ((float)p2->x[2]);

    posTemp[i*perTet + c++] = ((float)p4->x[0]);
    posTemp[i*perTet + c++] = ((float)p4->x[1]);
    posTemp[i*perTet + c++] = ((float)p4->x[2]);
  }
  return posTemp.data();
}

float* ParticleSystem::GetStrainAllTriColors(int* size, double strainSize) {
  *size = tets.size()* 4 * 3 * 3;
  int perTet = 4 * 3 * 3;
  colorTemp.resize(*size);
  for (int i = 0; i < tets.size(); i++) {
    Particle *p1,*p2,*p3,*p4;
    GetTetP(i, p1, p2, p3, p4);

    Eigen::Matrix3d temp;
    temp << p2->x - p1->x, p3->x - p1->x, p4->x - p1->x;
    Eigen::Matrix3d deformGradient = (temp * tets[i].inversePos) - Eigen::Matrix3d::Identity();
    Eigen::Matrix3d greenStrainTensor;
    if (corotational) {
      greenStrainTensor = .5 * (deformGradient + deformGradient.transpose() +
                                                deformGradient.transpose() * deformGradient);
    } else {
      greenStrainTensor = .5 * (deformGradient + deformGradient.transpose());
    }
    double v = .4;
    Eigen::VectorXd strainVec(6);
    strainVec << greenStrainTensor(0,0), greenStrainTensor(1,1),
                                    greenStrainTensor(2,2), greenStrainTensor(1,0),
                                    greenStrainTensor(1,2), greenStrainTensor(2,0);

    Eigen::MatrixXd strainToStress(6,6);
    strainToStress << 1 - v, v, v, 0, 0, 0,
                                           v, 1 - v, v, 0, 0, 0,
                                           v, v, 1 - v, 0, 0, 0,
                                           0, 0, 0, 1 - 2*v, 0, 0,
                                           0, 0, 0, 0, 1 - 2*v, 0,
                                           0, 0, 0, 0, 0, 1 - 2*v;
    Eigen::VectorXd stressVec(6);
    stressVec = (tets[i].k/((1 + v) * (1 - 2*v))) * strainToStress * strainVec;
    Eigen::Matrix3d stressTensor;
    stressTensor << stressVec[0], stressVec[3], stressVec[5],
                    stressVec[3], stressVec[1], stressVec[4],
                    stressVec[5], stressVec[4], stressVec[2];
    double strain = stressTensor.norm();
    tets[i].strain = strain;
  }
  for (int i = 0; i < tets.size(); i++) {
    double strain = tets[i].strain;
    int c = 0;
    while(c < perTet) {
      LerpColors(strain * strainSize, &(colorTemp[i*perTet+c]));
      c += 3;
    }
  }
  return colorTemp.data();
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

void ParticleSystem::SetupSingleSpring() {
  Reset();
  particles.emplace_back();
  particles.emplace_back();
  particles.emplace_back();
  //particles.emplace_back();
  particles[0].x << 0, 0, 0;
  particles[0].v << 0, 0.0, 0;
  particles[0].iMass = 1;
  particles[1].x << 1, .5, -1;
  particles[1].v << 0, 0.0, 0.0;
  particles[1].iMass = 1;
  particles[2].x << -1, .5, -1;
  particles[2].v << 0.0, 0.0, 0.0;
  particles[2].iMass = 1;

  //particles[3].x << 0, -.5, -1;
  //particles[3].v << 0.0, 0.0, 0.0;
  //particles[3].iMass = 1;

  CopyIntoStartPos();
  fixed_points.emplace_back();
  fixed_points[0].x << 0, -.5, -1;

  AddTet(0, 1, 2, -1);
  gravity = 2;
  ground = false;

}

void ParticleSystem::CopyIntoStartPos() {
  startPos.clear();
  for(int i = 0; i < particles.size(); ++i) {
    startPos.emplace_back();
    startPos[i] = particles[i].x;
  }
}

void ParticleSystem::SetupBendingBar() {
  Reset();
  int psize;
  double* points;
  std::vector<int> tets;
  MeshGen::GenerateBar(points, psize, tets, faces, facetotet);

  printf("Psize: %d, esize %d\n",psize, tets.size());

  for (int i = 0; i < psize; ++i) {
    particles.emplace_back();
    particles[i].x << points[i*3], points[i*3 + 1], points[i*3 + 2];
    particles[i].v << 0, 0, 0;
    particles[i].iMass = psize/20.0;
  }
  for (int i = 0; i < psize; ++i) {
    if (particles[i].x[2] == 0) {
      printf("fixed_point!\n");
      MakeFixedPoint(i, tets, faces);
      psize -= 1;
      i--;
    }
  }

  for (int i = 0; i < (tets.size()/4); ++i) {
    AddTet(tets[i*4], tets[i*4+1], tets[i*4 + 2], tets[i*4 + 3]);
  }
  //for (int i = 0; i < particles.size(); ++i) {
  //  CalculateParticleMass(i, 200.0/edges.size());
  //}
  CopyIntoStartPos(); 
  for (int i = 0; i < particles.size(); ++i) {
    //if (particles[i].x[2] < -2)
    //particles[i].v[1] += 5;
  }
  
  for (int i = 0; i < particles.size(); ++i) {
   // if (particles[i].x[2] >-2)
   // particles[i].v[1] += -5;
  }
  gravity = 9.8;
  ground = false;
  delete[] points;
}

void ParticleSystem::SetupArmadillo() {
  Reset();
  int psize;
  double* points;
  std::vector<int> tets;
  MeshGen::GenerateMesh(points, psize, tets, faces, facetotet, "Armadillo_simple2.ply");

  printf("Psize: %d, esize %d\n",psize, tets.size());

  for (int i = 0; i < psize; ++i) {
    particles.emplace_back();
    particles[i].x << points[i*3], points[i*3 + 1], points[i*3 + 2];
    particles[i].v << 0, 0, 0;
    particles[i].iMass = psize/20.0;
  }
  for (int i = 0; i < psize; ++i) {
    if (particles[i].x[1] < -6 && particles[i].x[0] < 2) {
      printf("fixed_point!\n");
      MakeFixedPoint(i, tets, faces);
      psize -= 1;
      i--;
    }
  }

  for (int i = 0; i < (tets.size()/4); ++i) {
    AddTet(tets[i*4], tets[i*4+1], tets[i*4 + 2], tets[i*4 + 3]);
  }
  //for (int i = 0; i < particles.size(); ++i) {
  //  CalculateParticleMass(i, 200.0/edges.size());
  //}
  CopyIntoStartPos(); 
  //for (int i = 0; i < particles.size(); ++i) {
  //  if (particles[i].x[2] < -6)
  //  particles[i].v[1] += .5;
  //}
  
  for (int i = 0; i < particles.size(); ++i) {
   // if (particles[i].x[2] >-2)
   // particles[i].v[1] += -5;
  }
  gravity = 9.8;
  //ground = true;
  delete[] points;
  printf("Number of faces%d\n", faces.size()/3);
}

void ParticleSystem::SetupMeshFile(char* filename) {
  Reset();
  int psize;
  double* points;
  std::vector<int> tets;
  MeshGen::GenerateMesh(points, psize, tets, faces, facetotet, filename);

  printf("Psize: %d, Number of tets %d\n",psize, tets.size());

  double lowestpoint = -10;
  for (int i = 0; i < psize; ++i) {
    particles.emplace_back();
    particles[i].x << points[i*3], points[i*3 + 1], points[i*3 + 2];
    particles[i].v << 0, 0, 0;
    particles[i].iMass = psize/20.0;
    if (points[i*3 + 1] > lowestpoint) {
      lowestpoint = points[i*3+1];
    }
  }

  //Make lowest points fixed
  for (int i = 0; i < psize; ++i) {
    if (particles[i].x[1] > lowestpoint - .1) {
      printf("fixed_point!\n");
      MakeFixedPoint(i, tets, faces);
      psize -= 1;
      i--;
    }
  }

  for (int i = 0; i < (tets.size()/4); ++i) {
    AddTet(tets[i*4], tets[i*4+1], tets[i*4 + 2], tets[i*4 + 3]);
  }
  CopyIntoStartPos(); 
  gravity = 9.8;
  ground = false;
  delete[] points;
  printf("Number of faces%d\n", faces.size()/3);
}

void ParticleSystem::MakeFixedPoint(int p, std::vector<int>& edges, std::vector<int>& faces) {
  fixed_points.emplace_back();
  fixed_points[fixed_points.size() - 1].x = particles[p].x;
  for (int i = 0; i < edges.size(); ++i) {
    if (edges[i] == p) {
      edges[i] = -1*fixed_points.size();
    } else if (edges[i] > p) {
      edges[i] -= 1;
    }
  }
  for (int i = 0; i < faces.size(); ++i) {
    if (faces[i] == p) {
      faces[i] = -1*fixed_points.size();
    } else if (faces[i] > p) {
      faces[i] -= 1;
    }
  }
  particles.erase(particles.begin() + p);
}

void ParticleSystem::CalculateParticleMass(int i, float springMass) {
 /* float mass = 0;
  for (int j = 0; j < springs.size(); ++j) {
    if (springs[j].to == i) mass += springMass/2;
    if (springs[j].from == i) mass += springMass/2;
  }
  if (mass == 0) mass = 1;
  particles[i].iMass = 1/mass;*/
}


void ParticleSystem::SetSpringProperties(double k, double c) {
  stiffness = k;
  dampness = c;
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
double equationSetupTime = 0;
double solveTime = 0;
};

void ParticleSystem::ImplicitEulerSparse(double timestep) {
  double tempTime = glfwGetTime();
  double curTime = tempTime;

  int vSize = 3 * particles.size();
  static Eigen::SparseMatrix<double> iesA;
  static Eigen::VectorXd iesb;

  static Eigen::SparseMatrix<double> iesdfdx;

  static std::vector<Eigen::Triplet<double>> iesdfdxtriplet;

  static std::vector<Eigen::Matrix3d> strainForTets;

  if (!hasPrev) {
    strainForTets.resize(tets.size() * 16);
    printf("Number of tets: %i\n", tets.size());
    for (int i = 0; i < tets.size(); i++) {
      Particle *p1,*p2,*p3,*p4;
      GetTetP(i, p1, p2, p3, p4);

      Eigen::Vector3d y0, y1, y2, y3;

      y1 << tets[i].inversePos(0,0), tets[i].inversePos(0,1), tets[i].inversePos(0,2);
      y2 << tets[i].inversePos(1,0), tets[i].inversePos(1,1), tets[i].inversePos(1,2);
      y3 << tets[i].inversePos(2,0), tets[i].inversePos(2,1), tets[i].inversePos(2,2);

      y0 = -1* y1 - y2 - y3;

      double v = .4;
      double a =  tets[i].posDet * tets[i].k * (1 - v) / ((1 + v) * (1 - 2 * v));
      double b =  tets[i].posDet * tets[i].k *  v / ((1 + v) * (1 - 2 * v));
      double c =  tets[i].posDet * tets[i].k * (1 - 2 * v) / ((1 + v) * (1 - 2 * v));
      Eigen::Matrix3d middle1, middle2;
      Eigen::Matrix3d temp, temp1,temp2,temp3,temp4;
      middle1 << a, b, b,
                 b, a, b,
                 b, b, a;
      middle2 << c, 0, 0,
                 0, c, 0,
                 0, 0, c;

      for (int index1 = 0; index1 < 4; ++index1) {
        Eigen::Vector3d *j0;
        switch(index1) {
          case 0: j0 = &y0; break;
          case 1: j0 = &y1; break;
          case 2: j0 = &y2; break;
          case 3: j0 = &y3; break;
        }
        temp1 << (*j0)[0], 0, 0,
                 0, (*j0)[1], 0,
                 0, 0, (*j0)[2];
        temp3 << (*j0)[1], 0, (*j0)[2],
                 (*j0)[0], (*j0)[2], 0,
                 0, (*j0)[1], (*j0)[0];
        for (int index2 = 0; index2 < 4; ++index2) {
          Eigen::Vector3d *j1;
          switch(index2) {
            case 0: j1 = &y0; break;
            case 1: j1 = &y1; break;
            case 2: j1 = &y2; break;
            case 3: j1 = &y3; break;
          }
          temp2 << (*j1)[0], 0, 0,
                   0, (*j1)[1], 0,
                   0, 0, (*j1)[2];
          temp4 << (*j1)[1], (*j1)[0], 0,
                   0, (*j1)[2], (*j1)[1],
                   (*j1)[2], 0, (*j1)[0];
          temp = temp1 * middle1 * temp2 + temp3 * middle2 * temp4;
          strainForTets[i * 16 + index1 * 4 + index2] = temp;
        }
      }
    }
  }

  iesA.resize(vSize, vSize);
  iesb.resize(vSize);
  iesdfdx.resize(vSize, vSize);

  iesdfdxtriplet.clear();

  Eigen::VectorXd f_0(vSize);
  f_0.setZero();

  for (int i = 0; i < tets.size(); i++) {
    Particle *p1,*p2,*p3,*p4;
    GetTetP(i, p1, p2, p3, p4);

    Eigen::Matrix3d Rot;
    if (corotational) {
      Eigen::Matrix3d m1,m2;
      Eigen::Vector3d r0,r1,r2;
      m1 << p2->x - p1->x, p3->x - p1->x, p4->x - p1->x;
      m2 = m1 * tets[i].inversePos;
      r0 = (m2.col(0)).normalized();
      r1 = (m2.col(1) - r0.dot(m2.col(1)) * r0).normalized();
      r2 = r0.cross(r1);
      Rot.col(0) = r0;
      Rot.col(1) = r1;
      Rot.col(2) = r2;

      //Eigen::Matrix3d mapping1, mapping2;
      //mapping1 << p2->x - p1->x, p3->x - p1->x, p4->x - p1->x;
      //mapping2 = mapping1 * tets[i].inversePos;
      //Eigen::Affine3d trans1;
      //trans1 = mapping2;
      //Rot = trans1.rotation();
    }
    // for all combos
    for (int index1 = 0; index1 < 4; ++index1) {
      for (int index2 = 0; index2 < 4; ++index2) {
        Eigen::Matrix3d* temp = &(strainForTets[i * 16 + index1 * 4 + index2]);
        Eigen::Matrix3d kelement;
        if (corotational) {
          if (tets[i].to[index1] >= 0) {
            kelement = Rot * (*temp) * Rot.transpose();
            Eigen::Vector3d oPos;


            if (tets[i].to[index2] >= 0) {
              oPos = startPos[tets[i].to[index2]];
              Eigen::Vector3d force = (Rot * (*temp) * oPos);
              f_0[tets[i].to[index1] * 3] += force[0];
              f_0[tets[i].to[index1] * 3 + 1] += force[1];
              f_0[tets[i].to[index1] * 3 + 2] += force[2];
              PushbackMatrix3d(iesdfdxtriplet, kelement, tets[i].to[index1] * 3, tets[i].to[index2] * 3, 1);
            } else {
              switch(index2) {
                case 0: oPos = p1->x; break;
                case 1: oPos = p2->x; break;
                case 2: oPos = p3->x; break;
                case 3: oPos = p4->x; break;
              }
              Eigen::Vector3d force = (Rot * (*temp) * oPos);
              force -= kelement * oPos; // since the fixed point doesn't exist in the solver
              f_0[tets[i].to[index1] * 3] += force[0];
              f_0[tets[i].to[index1] * 3 + 1] += force[1];
              f_0[tets[i].to[index1] * 3 + 2] += force[2];
            }
          }
        } else {
          if (tets[i].to[index1] < 0 || tets[i].to[index2] < 0) continue;
          kelement = *temp;
          PushbackMatrix3d(iesdfdxtriplet, kelement, tets[i].to[index1] * 3, tets[i].to[index2] * 3, 1);
        }
      }
    }
  }
  tempTime = glfwGetTime();
  tripletTime += tempTime - curTime;
  curTime = tempTime;

  iesdfdx.setFromTriplets(iesdfdxtriplet.begin(), iesdfdxtriplet.end());

  tempTime = glfwGetTime();
  fromTripletTime += tempTime - curTime;
  curTime = tempTime;

  Eigen::VectorXd v_0(vSize);
  Eigen::VectorXd x_0(vSize);
  Eigen::VectorXd f_ext(vSize);

  std::vector<Eigen::Triplet<double>> masstriplet;

  for (int i = 0; i < particles.size(); i++) {
    v_0[i * 3] = particles[i].v[0];
    v_0[i * 3 + 1] = particles[i].v[1];
    v_0[i * 3 + 2] = particles[i].v[2];
    if (corotational) {
      x_0[i * 3] = particles[i].x[0];
      x_0[i * 3 + 1] = particles[i].x[1];
      x_0[i * 3 + 2] = particles[i].x[2];
    } else {
      x_0[i * 3] = particles[i].x[0] - startPos[i][0];
      x_0[i * 3 + 1] = particles[i].x[1] - startPos[i][1];
      x_0[i * 3 + 2] = particles[i].x[2] - startPos[i][2];
    }
    f_ext[i * 3] = 0;
    f_ext[i * 3 + 1] = gravity/particles[i].iMass;
    f_ext[i * 3 + 2] = 0;
    masstriplet.push_back(Eigen::Triplet<double>(i*3,i*3,1/particles[i].iMass));
    masstriplet.push_back(Eigen::Triplet<double>(i*3+1,i*3+1,1/particles[i].iMass));
    masstriplet.push_back(Eigen::Triplet<double>(i*3+2,i*3+2,1/particles[i].iMass));
  }
  Eigen::VectorXd newv(vSize);
  //newv = v_0 + timestep * iesdfdx * x_0;

  iesA.setFromTriplets(masstriplet.begin(), masstriplet.end());

  iesb = iesA * v_0 + timestep * (iesdfdx * x_0 - f_0 + f_ext);
  iesA = iesA - (timestep * timestep * iesdfdx);
  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
  cg.setTolerance(.000001);
  cg.setMaxIterations(20);


  tempTime = glfwGetTime();
  equationSetupTime += tempTime - curTime;
  curTime = tempTime;

  cg.compute(iesA);
  if (hasPrev) newv = cg.solveWithGuess(iesb, vdiffprev);
  else newv = cg.solve(iesb);

  tempTime = glfwGetTime();
  solveTime += tempTime - curTime;
  curTime = tempTime;

  vdiffprev = newv;
  hasPrev = true;

  for (int i = 0; i < particles.size(); i++) {
    particles[i].v[0] = newv[i * 3];
    particles[i].v[1] = newv[i * 3 + 1];
    particles[i].v[2] = newv[i * 3 + 2];
    particles[i].x += timestep * particles[i].v;
  }
}

//void ParticleSystem::ImplicitEulerSparse(double timestep) {
//  int vSize = 3 * particles.size();
//  static Eigen::SparseMatrix<double> iesA;
//  static Eigen::VectorXd iesb;
//
//  static Eigen::SparseMatrix<double> iesdfdx;
//
//  static std::vector<Eigen::Triplet<double>> iesdfdxtriplet;
//
//  iesA.resize(vSize, vSize);
//  iesb.resize(vSize);
//
//  iesdfdxtriplet.clear();
//
//  if (!hasPrev) {
//    iesdfdx.resize(vSize, vSize);
//    for (int i = 0; i < tets.size(); i++) {
//      Particle *p1,*p2,*p3,*p4;
//      GetTetP(i, p1, p2, p3, p4);
//
//      Eigen::Vector3d y0, y1, y2, y3;
//
//      y1 << tets[i].inversePos(0,0), tets[i].inversePos(0,1), tets[i].inversePos(0,2);
//      y2 << tets[i].inversePos(1,0), tets[i].inversePos(1,1), tets[i].inversePos(1,2);
//      y3 << tets[i].inversePos(2,0), tets[i].inversePos(2,1), tets[i].inversePos(2,2);
//
//      y0 = -1* y1 - y2 - y3;
//
//      double v = .4;
//      double a =  tets[i].posDet * tets[i].k * (1 - v) / ((1 + v) * (1 - 2 * v));
//      double b =  tets[i].posDet * tets[i].k *  v / ((1 + v) * (1 - 2 * v));
//      double c =  tets[i].posDet * tets[i].k * (1 - 2 * v) / ((1 + v) * (1 - 2 * v));
//      Eigen::Matrix3d middle1, middle2;
//      Eigen::Matrix3d temp, temp1,temp2,temp3,temp4;
//      middle1 << a, b, b,
//                 b, a, b,
//                 b, b, a;
//      middle2 << c, 0, 0,
//                 0, c, 0,
//                 0, 0, c;
//
//      // for all combos
//      for (int index1 = 0; index1 < 4; ++index1) {
//        if (tets[i].to[index1] < 0) continue;
//        Eigen::Vector3d *j0;
//        switch(index1) {
//          case 0: j0 = &y0; break;
//          case 1: j0 = &y1; break;
//          case 2: j0 = &y2; break;
//          case 3: j0 = &y3; break;
//        }
//        temp1 << (*j0)[0], 0, 0,
//                 0, (*j0)[1], 0,
//                 0, 0, (*j0)[2];
//        temp3 << (*j0)[1], 0, (*j0)[2],
//                 (*j0)[0], (*j0)[2], 0,
//                 0, (*j0)[1], (*j0)[0];
//        for (int index2 = 0; index2 < 4; ++index2) {
//          if (tets[i].to[index2] >= 0) {//continue;
//          Eigen::Vector3d *j1;
//          switch(index2) {
//            case 0: j1 = &y0; break;
//            case 1: j1 = &y1; break;
//            case 2: j1 = &y2; break;
//            case 3: j1 = &y3; break;
//          }
//          temp2 << (*j1)[0], 0, 0,
//                   0, (*j1)[1], 0,
//                   0, 0, (*j1)[2];
//          temp4 << (*j1)[1], (*j1)[0], 0,
//                   0, (*j1)[2], (*j1)[1],
//                   (*j1)[2], 0, (*j1)[0];
//          temp = temp1 * middle1 * temp2 + temp3 * middle2 * temp4;
//          PushbackMatrix3d(iesdfdxtriplet, temp, tets[i].to[index1] * 3, tets[i].to[index2] * 3, 1);
//        } else {}}
//      }
//    }
//    iesdfdx.setFromTriplets(iesdfdxtriplet.begin(), iesdfdxtriplet.end());
//  }
//
//  Eigen::VectorXd v_0(vSize);
//  Eigen::VectorXd x_0(vSize);
//  Eigen::VectorXd f_ext(vSize);
//
//  std::vector<Eigen::Triplet<double>> masstriplet;
//
//  for (int i = 0; i < particles.size(); i++) {
//    v_0[i * 3] = particles[i].v[0];
//    v_0[i * 3 + 1] = particles[i].v[1];
//    v_0[i * 3 + 2] = particles[i].v[2];
//    x_0[i * 3] = particles[i].x[0] - startPos[i][0];
//    x_0[i * 3 + 1] = particles[i].x[1] - startPos[i][1];
//    x_0[i * 3 + 2] = particles[i].x[2] - startPos[i][2];
//    f_ext[i * 3] = 0;
//    f_ext[i * 3 + 1] = gravity/particles[i].iMass;
//    f_ext[i * 3 + 2] = 0;
//    masstriplet.push_back(Eigen::Triplet<double>(i*3,i*3,1/particles[i].iMass));
//    masstriplet.push_back(Eigen::Triplet<double>(i*3+1,i*3+1,1/particles[i].iMass));
//    masstriplet.push_back(Eigen::Triplet<double>(i*3+2,i*3+2,1/particles[i].iMass));
//  }
//  Eigen::VectorXd newv(vSize);
//  //newv = v_0 + timestep * iesdfdx * x_0;
//
//  iesA.setFromTriplets(masstriplet.begin(), masstriplet.end());
//  iesb = iesA * v_0 + timestep * (iesdfdx * x_0 + f_ext);
//  iesA = iesA - (timestep * timestep * iesdfdx);
//  Eigen::ConjugateGradient<Eigen::SparseMatrix<double> > cg;
//  cg.setTolerance(.001);
//  cg.setMaxIterations(30);
//
//  double tempTime = glfwGetTime();
//  double curTime = tempTime;
//
//  cg.compute(iesA);
//  if (hasPrev) newv = cg.solveWithGuess(iesb, vdiffprev);
//  else newv = cg.solve(iesb);
//
//  tempTime = glfwGetTime();
//  solveTime += tempTime - curTime;
//  curTime = tempTime;
//
//  vdiffprev = newv;
//  hasPrev = true;
//
//  for (int i = 0; i < particles.size(); i++) {
//    particles[i].v[0] = newv[i * 3];
//    particles[i].v[1] = newv[i * 3 + 1];
//    particles[i].v[2] = newv[i * 3 + 2];
//    particles[i].x += timestep * particles[i].v;
//  }
//}

/*void ParticleSystem::ComputeForces() {
  //Zero all forces
  for (int i = 0; i < particles.size(); i++) {
    particles[i].f << 0.0, gravity/particles[i].iMass, 0.0;
  }

  //Compute spring forces
  int sSize = tets.size();
  for (int i = 0; i < sSize; i++) {
    Particle *p1,*p2,*p3,*p4;
    GetTetP(i, p1, p2, p3, p4);

    Eigen::Matrix3d temp;
    temp << p2->x - p1->x, p3->x - p1->x, p4->x - p1->x;
    Eigen::Matrix3d deformGradient = (temp * tets[i].inversePos) - Eigen::Matrix3d::Identity();


    Eigen::Matrix3d greenStrainTensor = .5 * (deformGradient + deformGradient.transpose() +
                                              deformGradient.transpose() * deformGradient);
    std::cout << greenStrainTensor << std::endl;
    double v = .4;
    Eigen::VectorXd strainVec(6);
    strainVec << greenStrainTensor(0,0), greenStrainTensor(1,1),
                                    greenStrainTensor(2,2), greenStrainTensor(1,0),
                                    greenStrainTensor(1,2), greenStrainTensor(2,0);

    Eigen::MatrixXd strainToStress(6,6);
    strainToStress << 1 - v, v, v, 0, 0, 0,
                                           v, 1 - v, v, 0, 0, 0,
                                           v, v, 1 - v, 0, 0, 0,
                                           0, 0, 0, 1 - 2*v, 0, 0,
                                           0, 0, 0, 0, 1 - 2*v, 0,
                                           0, 0, 0, 0, 0, 1 - 2*v;
    Eigen::VectorXd stressVec(6);
    stressVec = (tets[i].k/((1 + v) * (1 - 2*v))) * strainToStress * strainVec;
    Eigen::Matrix3d stressTensor;
    stressTensor << stressVec[0], stressVec[3], stressVec[5],
                    stressVec[3], stressVec[1], stressVec[4],
                    stressVec[5], stressVec[4], stressVec[2];
    std::cout << "STRESS" << std::endl;
    std::cout << stressTensor << std::endl;
    tets[i].strain = stressTensor.norm();
    //stressTensor = greenStrainTensor;
    // Lets loop over faces and apply stress
    for(int j = 0; j < 4; j++) {
      Particle *j0, *j1, *j2;
      switch(j) {
        case 0:
          j0 = p1; j1 = p3; j2 = p2; break;
        case 1:
          j0 = p1; j1 = p2; j2 = p4; break;
        case 2:
          j0 = p1; j1 = p4; j2 = p3; break;
        case 3:
          j0 = p3; j1 = p4; j2 = p2; break;
       }
       Eigen::Vector3d faceForce = stressTensor * (j1->x - j0->x).cross(j2->x - j0->x);
       printf("face force %d , %f\n", j, faceForce.norm());
       printf("COMface force com %f %f %f\n", j, faceForce[0], faceForce[1], faceForce[2]);
       j0->f += 1.0/3.0 * faceForce;
       j1->f += 1.0/3.0 * faceForce;
       j2->f += 1.0/3.0 * faceForce;
    }
  }
  printf("p1 f %f\n",particles[0].f[0]);
}*/
void ParticleSystem::ComputeForces() {}
void ParticleSystem::ExplicitEuler(double timestep) {
  /*phaseTemp.resize(particles.size() * 6);
  for (int i = 0; i < particles.size(); i++) {
    phaseTemp[i * 6] = particles[i].x[0];
    phaseTemp[i * 6 + 1] = particles[i].x[1];
    phaseTemp[i * 6 + 2] = particles[i].x[2];
    phaseTemp[i * 6 + 3] = particles[i].v[0];
    phaseTemp[i * 6 + 4] = particles[i].v[1];
    phaseTemp[i * 6 + 5] = particles[i].v[2];
  }*/
  ComputeForces();
  for (int i = 0; i < particles.size(); i++) {
    particles[i].v += particles[i].f * particles[i].iMass * timestep;
    particles[i].x += particles[i].v * timestep;
  }
  /*
  ComputeForces();
  for (int i = 0; i < particles.size(); i++) {
    particles[i].x[0] = phaseTemp[i * 6];
    particles[i].x[1] = phaseTemp[i * 6 + 1];
    particles[i].x[2] = phaseTemp[i * 6 + 2];
    particles[i].v[0] = phaseTemp[i * 6 + 3];
    particles[i].v[1] = phaseTemp[i * 6 + 4];
    particles[i].v[2] = phaseTemp[i * 6 + 5];
    particles[i].v += particles[i].f * particles[i].iMass * timestep;
    //particles[i].v *= .95;
    particles[i].x += particles[i].v * timestep;
  }*/
}

void ParticleSystem::GetProfileInfo(double& triplet, double& fromTriplet, double& solve, double& setupTime) {
   triplet = tripletTime;
   fromTriplet = fromTripletTime;
   solve = solveTime;
   setupTime = equationSetupTime;
}

void ParticleSystem::AddTet(int x1, int x2, int x3, int x4) {
  tets.emplace_back();
  int index = tets.size() - 1;
  tets[index].to[0] = x1;
  tets[index].to[1] = x2;
  tets[index].to[2] = x3;
  tets[index].to[3] = x4;

  tets[index].k = stiffness;
  tets[index].c = dampness;

  Particle *p1,*p2,*p3,*p4;
  GetTetP(index, p1, p2, p3, p4);
  Eigen::Matrix3d temp;
  temp << p2->x - p1->x, p3->x - p1->x, p4->x - p1->x;

  tets[index].oldPos[0] = p2->x - p1->x;
  tets[index].oldPos[1] = p3->x - p1->x;
  tets[index].oldPos[2] = p4->x - p1->x;
  tets[index].posDet = temp.determinant();
  tets[index].inversePos = temp.inverse();
}

void ParticleSystem::GetTetP(int i, Particle*& x1, Particle*& x2, Particle*& x3, Particle*& x4) {
  if (tets[i].to[0] < 0)
    x1 = &(fixed_points[tets[i].to[0] * -1 - 1]);
  else
    x1 = &(particles[tets[i].to[0]]);

  if (tets[i].to[1] < 0)
    x2 = &(fixed_points[tets[i].to[1] * -1 - 1]);
  else
    x2 = &(particles[tets[i].to[1]]);

  if (tets[i].to[2] < 0)
    x3 = &(fixed_points[tets[i].to[2] * -1 - 1]);
  else
    x3 = &(particles[tets[i].to[2]]);

  if (tets[i].to[3] < 0)
    x4 = &(fixed_points[tets[i].to[3] * -1 - 1]);
  else
    x4 = &(particles[tets[i].to[3]]);
}

void ParticleSystem::GetPointP(int i, Particle*& x1) {
  if (i < 0)
    x1 = &(fixed_points[i * -1 - 1]);
  else
    x1 = &(particles[i]);
}

void ParticleSystem::Reset() {
  particles.clear();
  fixed_points.clear();
  tets.clear();
  hasPrev = false;
  faces.clear();
  facetotet.clear();
}
