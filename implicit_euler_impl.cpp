#include <GLFW/glfw3.h>
#include "particle_system.h"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"

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

void ParticleSystem::GetProfileInfo(double& triplet, double& fromTriplet, double& solve, double& setupTime) {
   triplet = tripletTime;
   fromTriplet = fromTripletTime;
   solve = solveTime;
   setupTime = equationSetupTime;
}

void ParticleSystem::Reset() {
  particles.clear();
  fixed_points.clear();
  tets.clear();
  hasPrev = false;
  faces.clear();
  facetotet.clear();
}

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

      double v = volConserve;
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
    f_ext[i * 3] = particles[i].f[0];
    f_ext[i * 3 + 1] = gravity/particles[i].iMass + particles[i].f[1];
    f_ext[i * 3 + 2] = particles[i].f[2];
    particles[i].f << 0,0,0;
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
