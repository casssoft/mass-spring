#include "particle_system.h"
#include "Eigen/Sparse"
#include "Eigen/Dense"
#include "Eigen/IterativeLinearSolvers"
#include <iostream>
#include <math.h>

#ifdef COLLISION_SELFCCD
#include "collision_system.h"
#endif

#ifdef COLLISION_PQP
#include "collision_system_pqp.h"
static int initialFaceSize;
#endif
static std::vector<int> faceToOut;
void ParticleSystem::HandleCollisions(double timestep) {
  faceToOut.clear();
  if (useColSys) {
#ifdef COLLISION_SELFCCD
    for (int i = 0; i < outsidePoints.size(); i++) {
      if (outsidePoints[i] >= 0) {
        colSys->UpdateVertex(i, particles[outsidePoints[i]].x);
      }
    }
    std::vector<unsigned int> vertexToFace;
    std::vector<unsigned int> edgeToEdge;
    std::vector<float> veToFaTime;
    std::vector<float> edToEdTime;
    colSys->GetCollisions(vertexToFace, edgeToEdge, veToFaTime, edToEdTime);
    if (colRolBack) {
      int colCount = 0;
      int earliestIndex = 0;
      double curTimeStep = timestep;
      while(earliestIndex >= 0) {
        //printf("Start loop\n");
        earliestIndex = -1;
        double eTime = 0;
        for (int i = 0; i < vertexToFace.size(); i += 2) {
          // calculate normal of tri
          int p1_i, p2_i, p3_i, v1_i;
          p1_i = outsidePoints[faceToOut[3 * vertexToFace[i + 1]]];
          p2_i = outsidePoints[faceToOut[3 * vertexToFace[i + 1] + 1]];
          p3_i = outsidePoints[faceToOut[3 * vertexToFace[i + 1] + 2]];
          v1_i = outsidePoints[vertexToFace[i]];
          if (p1_i < 0 && v1_i >= 0) {
            if (earliestIndex == -1 || veToFaTime[i/2] < eTime) {
              eTime = veToFaTime[i/2];
              earliestIndex = i;
              //printf("Found early collision\n");
            }
          }
        }
        if (earliestIndex >= 0) {
            Particle *p1, *p2, *p3, *v1;
            int p1_i, p2_i, p3_i, v1_i;
            p1_i = outsidePoints[faceToOut[3 * vertexToFace[earliestIndex + 1]]];
            p2_i = outsidePoints[faceToOut[3 * vertexToFace[earliestIndex + 1] + 1]];
            p3_i = outsidePoints[faceToOut[3 * vertexToFace[earliestIndex + 1] + 2]];
            v1_i = outsidePoints[vertexToFace[earliestIndex]];
            GetPointP(p1_i, p1);
            GetPointP(p2_i, p2);
            GetPointP(p3_i, p3);
            GetPointP(v1_i, v1);
            Eigen::Vector3d temp1, temp2;
            temp1 = p2->x - p1->x;
            temp2 = p3->x - p1->x;
            temp1 = temp1.cross(temp2);
            temp1.normalize();
            temp1;
            // Project vertex onto plane
            double d = p1->x.dot(temp1);
            double v = (d - (v1->x.dot(temp1)));
            if (v < 0) {
              //printf("inside\n");
            }
            Eigen::Vector3d planePoint = v1->x + v * temp1;
            //printf("Original point: %f, %f, %f\n", v1->x[0], v1->x[1], v1->x[2]);
            //printf("Plane point: %f, %f, %f\n", planePoint[0], planePoint[1], planePoint[2]);
            //printf("V: %f\n", v);
            //printf("temp1: %f, %f, %f\n", temp1[0], temp1[1], temp1[2]);
            planePoint +=  temp1 * .05 * 60 * curTimeStep * (1 + .0001 * colCount);

            for (int i = 0; i < particles.size(); i++) {
              Eigen::Vector3d temp;
              temp = eTime * particles[i].x + (1 - eTime) * prevPos[i];
              particles[i].x = temp;
              temp = eTime * particles[i].v + (1 - eTime) * prevVel[i];
              particles[i].v = temp;
              temp = eTime * particles[i].f + (1 - eTime) * prevFEXT[i];
              particles[i].f = temp;
            }
            v1->x[0] = planePoint[0];
            v1->x[1] = planePoint[1];
            v1->x[2] = planePoint[2];
            v1->v[0] = 0;
            v1->v[1] = 0;
            v1->v[2] = 0;

            for (int i = 0; i < outsidePoints.size(); i++) {
              if (outsidePoints[i] >= 0) {
                colSys->UpdateVertex(i, particles[outsidePoints[i]].x);
              }
            }
            colSys->GetCollisions(vertexToFace, edgeToEdge, veToFaTime, edToEdTime);
            vertexToFace.clear();
            edgeToEdge.clear();
            veToFaTime.clear();
            edToEdTime.clear();
            //printf("Doing another round of euler with time step %f\n", curTimeStep * (1 - eTime));
            ImplicitEulerSparse(curTimeStep * (1 - eTime));
            curTimeStep -= eTime * curTimeStep;
            for (int i = 0; i < outsidePoints.size(); i++) {
              if (outsidePoints[i] >= 0) {
                colSys->UpdateVertex(i, particles[outsidePoints[i]].x);
              }
            }
            colSys->GetCollisions(vertexToFace, edgeToEdge, veToFaTime, edToEdTime);
            colCount += 1;
            //fprintf(stderr, "c: %i curTimeStep: %f ", colCount, curTimeStep);
            if (colCount > 60) break;
        }
      }
      fprintf(stderr, "ColCount: %i\n", colCount);
    } else {
    for (int i = 0; i < vertexToFace.size(); i += 2) {
      // calculate normal of tri
      Particle *p1, *p2, *p3, *v1;
      int p1_i, p2_i, p3_i, v1_i;
      p1_i = outsidePoints[faceToOut[3 * vertexToFace[i + 1]]];
      p2_i = outsidePoints[faceToOut[3 * vertexToFace[i + 1] + 1]];
      p3_i = outsidePoints[faceToOut[3 * vertexToFace[i + 1] + 2]];
      v1_i = outsidePoints[vertexToFace[i]];
      if (p1_i < 0 && v1_i >= 0) {
        GetPointP(p1_i, p1);
        GetPointP(p2_i, p2);
        GetPointP(p3_i, p3);
        GetPointP(v1_i, v1);
        Eigen::Vector3d temp1, temp2;
        temp1 = p2->x - p1->x;
        temp2 = p3->x - p1->x;
        temp1 = temp1.cross(temp2);
        temp1.normalize();
        temp1;
        // Project vertex onto plane
        double d = p1->x.dot(temp1);
        double v = (d - (v1->x.dot(temp1)));
        if (v < 0) {
          //printf("inside\n");
        }
        Eigen::Vector3d planePoint = v1->x + v * temp1;
        //printf("Original point: %f, %f, %f\n", v1->x[0], v1->x[1], v1->x[2]);
        //printf("Plane point: %f, %f, %f\n", planePoint[0], planePoint[1], planePoint[2]);
        //printf("V: %f\n", v);
        //printf("temp1: %f, %f, %f\n", temp1[0], temp1[1], temp1[2]);
        planePoint +=  temp1 * .05 * 30 * timestep * timestep;
        v1->x[0] = planePoint[0];
        v1->x[1] = planePoint[1];
        v1->x[2] = planePoint[2];
        //v1->x[0] -= v1->v[0] * timestep;
        //v1->x[1] -= v1->v[1] * timestep;
        //v1->x[2] -= v1->v[2] * timestep;
        v1->v[0] = 0;
        v1->v[1] = 0;
        v1->v[2] = 0;
        colSys->UpdateVertex(vertexToFace[i], v1->x);
      }
    }
    for (int i = 0; i < edgeToEdge.size(); i += 4) {
      // calculate normal of tri
      Particle *p1, *p2, *p3, *p4;
      int p1_i, p2_i, p3_i, p4_i;
      p1_i = outsidePoints[edgeToEdge[i]];
      p2_i = outsidePoints[edgeToEdge[i + 1]];
      p3_i = outsidePoints[edgeToEdge[i + 2]];
      p4_i = outsidePoints[edgeToEdge[i + 3]];
      // is one edge fixed points and the other not?
      bool switched = false;
      if (p1_i < 0 && p2_i < 0 && p3_i >= 0 && p4_i >= 0) {
        int temp = p1_i;
        p1_i = p3_i;
        p3_i = temp;

        temp = p2_i;
        p2_i = p4_i;
        p4_i = temp;
        switched = true;
      }

      if (p1_i >= 0 && p2_i >= 0 && p3_i < 0 && p4_i < 0) {
        GetPointP(p1_i, p1);
        GetPointP(p2_i, p2);
        GetPointP(p3_i, p3);
        GetPointP(p4_i, p4);

        // Get col point and put it right before?
        //Eigen::Vector3d newPos;
        //float time = edToEdTime[i/4];
        //if (time < .2) time = 0;
        //else time -= .2;
        //time = 0;
        //newPos = p1->x * time + prevPos[p1_i] * (1 - time);
        //p1->x = newPos;
        //newPos = p2->x * time + prevPos[p2_i] * (1 - time);
        //p2->x = newPos;
        //p1->v << 0,0,0;
        //p2->v << 0,0,0;
        //if (switched) {
        //  colSys->UpdateVertex(edgeToEdge[i + 2], p1->x);
        //  colSys->UpdateVertex(edgeToEdge[i + 3], p2->x);
        //} else {
        //  colSys->UpdateVertex(edgeToEdge[i], p1->x);
        //  colSys->UpdateVertex(edgeToEdge[i + 1], p2->x);
        //}
      }
    }
    vertexToFace.clear();
    edgeToEdge.clear();
    veToFaTime.clear();
    edToEdTime.clear();
    colSys->GetCollisions(vertexToFace, edgeToEdge, veToFaTime, edToEdTime);

//Eigen::Vector3d newVelocity = (p1->v * bary[0] + p2->v * bary[1] + p3->v * bary[2]) / 2;
//v1->x = planePoint;
//v1->v = newVelocity;

      //if (p1_i >= 0 && p2_i >= 0 && p3_i >= 0 && v1_i >= 0) {
      //  p1 = &(particles[p1_i]);
      //  p2 = &(particles[p2_i]);
      //  p3 = &(particles[p3_i]);
      //  v1 = &(particles[v1_i]);
      //  Eigen::Vector3d temp1, temp2;
      //  temp1 = p2->x - p1->x;
      //  temp2 = p3->x - p1->x;
      //  temp1 = temp1.cross(temp2);
      //  temp1.normalize();
      //  temp1;
      //  // Project vertex onto plane
      //  double d = p1->x.dot(temp1);
      //  double v = (d - (v1->x.dot(temp1)));
      //  Eigen::Vector3d planePoint = v1->x + v * temp1;
      //  printf("Original point: %f, %f, %f\n", v1->x[0], v1->x[1], v1->x[2]);
      //  printf("Plane point: %f, %f, %f\n", planePoint[0], planePoint[1], planePoint[2]);
      //  printf("V: %f\n", v);
      //  printf("temp1: %f, %f, %f\n", temp1[0], temp1[1], temp1[2]);
      //  Eigen::Matrix3d baryM;
      //  baryM << p1->x[0], p2->x[0], p3->x[0],
      //           p1->x[1], p2->x[1], p3->x[1],
      //           p1->x[2], p2->x[2], p3->x[2];

      //  Eigen::Vector3d bary = baryM.inverse() * planePoint;
      //  printf("Sum of barys: %f\n", bary[0] + bary[1] + bary[2]);
      //  Eigen::Vector3d newVelocity = (v1->v + p1->v * bary[0] + p2->v * bary[1] + p3->v * bary[2]) / 2;
      //  v1->x = (v1->x + planePoint)/2;
      //  Eigen::Vector3d moveTri = 3 * (v1->x - planePoint);
      //  v1->v = newVelocity;
      //  p1->v = newVelocity * bary[0] + p1->v * (1 - bary[0]);
      //  p2->v = newVelocity * bary[1] + p2->v * (1 - bary[1]);
      //  p3->v = newVelocity * bary[2] + p3->v * (1 - bary[2]);
      //  p1->x += moveTri * bary[0];
      //  p2->x += moveTri * bary[1];
      //  p3->x += moveTri * bary[2];

      //  //Eigen::Vector3d newVelocity = (p1->v * bary[0] + p2->v * bary[1] + p3->v * bary[2]) / 2;
      //  //v1->x = planePoint;
      //  //v1->v = newVelocity;

      //  //colSys->UpdateVertex(vertexToFace[i], v1->x);
      //  //colSys->UpdateVertex(faceToOut[vertexToFace[i + 1]], p1->x);
      //  //colSys->UpdateVertex(faceToOut[vertexToFace[i + 1] + 1], p2->x);
      //  //colSys->UpdateVertex(faceToOut[vertexToFace[i + 1] + 2], p3->x);
      //  // for now just move the vertex to match the face
      //  // calculate barycentric coordinates of hit
      //}
  //  }
  //}
  }
  }
#endif // COLLISION_SELFCCD
#ifdef COLLISION_PQP
  // reinit object
  std::vector<Eigen::Vector3d> verts;
  std::vector<int> otris;
  for (int i = 0; i < initialFaceSize; ++i) {
    Particle*x;
    GetPointP(faces[i], x);
    verts.push_back(x->x);
    otris.push_back(i);
  }
  colSys->InitObjectModel(verts, otris);

  std::vector<unsigned int> vertexToFace;
  std::vector<unsigned int> edgeToEdge;
  std::vector<Eigen::Vector3d> moveEdge;
  std::vector<double> edgeU;
  colSys->GetCollisions(vertexToFace, edgeToEdge, edgeU, moveEdge);
  for (int i = 0; i < edgeToEdge.size(); i += 2) {
    break;
    Particle *v1, *v2;
    int v1_i, v2_i;
    v1_i = faces[edgeToEdge[i]];
    v2_i = faces[edgeToEdge[i + 1]];
    //fprintf(stderr, "got collisions with face %i %i %i and vertex %i\n", p1_i, p2_i, p3_i, v1_i);
    if (v1_i >= 0 && v2_i >= 0) {
      GetPointP(v1_i, v1);
      GetPointP(v2_i, v2);
      //fprintf(stderr, "v1_i %i v2_i %i prevPos size %i\n", v1_i, v2_i, prevPos.size());
      double u = edgeU[i/2];
      double mu = 1.0 / (fabs(.5 - u) + .5);
      v1->x = prevPos[v1_i];
      v2->x = prevPos[v2_i];
      Eigen::Matrix<double, 9, 9> m;
      m << 1, 0, 0,  0, 0, 0,  1 - u, 0, 0,
           0, 1, 0,  0, 0, 0,  0, 1 - u, 0,
           0, 1, 1,  0, 0, 0,  0, 0, 1 - u,

           0, 0, 0,  1, 0, 0,  u, 0, 0,
           0, 0, 0,  0, 1, 0,  0, u, 0,
           0, 0, 0,  0, 0, 1,  0, 0, u,

           1 - u, 0, 0,  u, 0, 0,  0, 0, 0,
           0, 1 - u, 0,  0, u, 0,  0, 0, 0,
           0, 0, 1 - u,  0, 0, u,  0, 0, 0;
      //fprintf(stderr, "assigned m\n");
      Eigen::VectorXd x(9), b(9);
      b << v1->v[0], v1->v[1], v1->v[2], v2->v[0], v2->v[1], v2->v[2], 0, 0, 0;
      //fprintf(stderr, "assigned b\n");
      x = m.colPivHouseholderQr().solve(b);
      //fprintf(stderr, "after householder\n", v1_i, v2_i, prevPos.size());
      v1->v[0] = x(0);
      v1->v[1] = x(1);
      v1->v[2] = x(2);
      v2->v[0] = x(3);
      v2->v[1] = x(4);
      v2->v[2] = x(5);
      v1->mark = true;
      v2->mark = true;
    }
  }
  for (int i = 0; i < vertexToFace.size(); i += 2) {
      Particle *p1, *p2, *p3, *v1;
      int p1_i, p2_i, p3_i, v1_i;
      p1_i = faces[vertexToFace[i + 1] + initialFaceSize];
      p2_i = faces[vertexToFace[i + 1] + 1 + initialFaceSize];
      p3_i = faces[vertexToFace[i + 1] + 2 + initialFaceSize];
      v1_i = faces[vertexToFace[i]];
      //fprintf(stderr, "got collisions with face %i %i %i and vertex %i\n", p1_i, p2_i, p3_i, v1_i);
      if (p1_i < 0 && v1_i >= 0) {
        GetPointP(p1_i, p1);
        GetPointP(p2_i, p2);
        GetPointP(p3_i, p3);
        GetPointP(v1_i, v1);
        Eigen::Vector3d temp1, temp2;
        temp1 = p2->x - p1->x;
        temp2 = p3->x - p1->x;
        temp1 = temp1.cross(temp2);
        if (temp1.norm() == 0) {
          continue;
        }
        temp1.normalize();
        // Project vertex onto plane
        double d = p1->x.dot(temp1);
        double v = (d - (v1->x.dot(temp1)));
        if (v < 0) {
          //printf("inside\n");
        }
        Eigen::Vector3d planePoint = v1->x + v * temp1;
        //printf("Original point: %f, %f, %f\n", v1->x[0], v1->x[1], v1->x[2]);
        //printf("Plane point: %f, %f, %f\n", planePoint[0], planePoint[1], planePoint[2]);
        //printf("V: %f\n", v);
        //printf("temp1: %f, %f, %f\n", temp1[0], temp1[1], temp1[2]);
        planePoint +=  temp1 * .05 * timestep;
        v1->x[0] = planePoint[0];
        v1->x[1] = planePoint[1];
        v1->x[2] = planePoint[2];
        //v1->x[0] -= v1->v[0] * timestep;
        //v1->x[1] -= v1->v[1] * timestep;
        //v1->x[2] -= v1->v[2] * timestep;
        v1->v[0] = 0;
        v1->v[1] = 0;
        v1->v[2] = 0;
        //v1->mark = true;
      }
    }
  }
#endif 
}
void ParticleSystem::SetupCollisions(double lowestpoint) {
  initialFaceSize = faces.size();
  Particle *g1,* g2,* g3, * g4, *g5, *g6, *g7, *g8;
  g1 = new Particle();
  g2 = new Particle();
  g3 = new Particle();
  g4 = new Particle();
  g5 = new Particle();
  g6 = new Particle();
  g7 = new Particle();
  g8 = new Particle();
  g1->x << 2, lowestpoint + 1, 2;
  g1->v << 0, 0, 0;
  g2->x << -2, lowestpoint + 1, 2;
  g2->v << 0, 0, 0;
  g3->x << -2, lowestpoint + 1, -2;
  g3->v << 0, 0, 0;
  g4->x << 2, lowestpoint + 1, -2;
  g4->v << 0, 0, 0;


  g5->x << 2, lowestpoint + 5, 2;
  g5->v << 0, 0, 0;
  g6->x << -2, lowestpoint + 5, 2;
  g6->v << 0, 0, 0;
  g7->x << -2, lowestpoint + 5, -2;
  g7->v << 0, 0, 0;
  g8->x << 2, lowestpoint + 5, -2;
  g8->v << 0, 0, 0;

  fixed_points.push_back(*g1);
  fixed_points.push_back(*g2);
  fixed_points.push_back(*g3);
  fixed_points.push_back(*g4);
  fixed_points.push_back(*g5);
  fixed_points.push_back(*g6);
  fixed_points.push_back(*g7);
  fixed_points.push_back(*g8);
  delete g1;
  delete g2;
  delete g3;
  delete g4;
  delete g5;
  delete g6;
  delete g7;
  delete g8;
  //outsidePoints.push_back(-1 * fixed_points.size()); // 1, 3, -1
  //outsidePoints.push_back(-1 * fixed_points.size() + 1); // -1, 3 -1
  //outsidePoints.push_back(-1 * fixed_points.size() + 2); // -1, 3, 1
  //outsidePoints.push_back(-1 * fixed_points.size() + 3); // 1, 3, 1
  //outsidePoints.push_back(-1 * fixed_points.size() + 4); // 1, 1, -1
  //outsidePoints.push_back(-1 * fixed_points.size() + 5); // -1, 1, -1
  //outsidePoints.push_back(-1 * fixed_points.size() + 6); // -1, 1, 1
  //outsidePoints.push_back(-1 * fixed_points.size() + 7); // 1, 1, 1

  //  6 ____7
  //  / \  / \ <--- 3
  //2/___ \5___\4
  //  \   /\   /
  //   1\/___\/0
  //
  //
#define F_PUSH(num) faces.push_back(-1 * fixed_points.size() + num);
#define F_SIDE(v1, v2, v3, v4) \
  F_PUSH(v1) F_PUSH(v2) F_PUSH(v4) \
  F_PUSH(v2) F_PUSH(v3) F_PUSH(v4)

  F_SIDE(7, 6, 5, 4) //top
  F_SIDE(4, 5, 1, 0)//front
  F_SIDE(0, 1, 2, 3) //bottom
  F_SIDE(6, 7, 3, 2) // back
  F_SIDE(7, 4, 0, 3) // right side
  F_SIDE(5, 6, 2, 1) // left side

  useColSys = true;
  outsidePoints.clear();

  faceToOut.clear();
  for (int i = 0; i < faces.size(); ++i) {
    faceToOut.push_back(0);
    bool addPoint = true;
    //printf("loop %d\n", i);
    for(int j = outsidePoints.size() - 1; j >= 0 && addPoint; --j) {
      if (outsidePoints[j] == faces[i]) {
        faceToOut[faceToOut.size() - 1] = j;
        addPoint = false;
      }
    }
    if (addPoint) {
      outsidePoints.push_back(faces[i]);
      faceToOut[faceToOut.size() - 1] = outsidePoints.size() - 1;
    }
  }

#ifdef COLLISION_SELFCCD
  std::vector<Eigen::Vector3d> verts;
  for (int i = 0; i < outsidePoints.size(); ++i) {
    Particle* x;
    //printf("Getting point %i\n", outsidePoints[i]);
    GetPointP(outsidePoints[i], x);
    verts.push_back(x->x);
  }
  //printf("starting init System\n");
  colSys->InitSystem(verts, faceToOut);
#endif
#ifdef COLLISION_PQP
  std::vector<Eigen::Vector3d> verts;
  std::vector<int> otris;
  for (int i = 0; i < initialFaceSize; ++i) {
    Particle*x;
    GetPointP(faces[i], x);
    verts.push_back(x->x);
    otris.push_back(i);
  }
  printf("size of verts %i, size of otris %i\n", verts.size(), otris.size());
  std::vector<Eigen::Vector3d> groundverts;
  std::vector<int> gtris;
  for (int i = initialFaceSize; i < faces.size(); ++i) {
    Particle*x;
    GetPointP(faces[i], x);
    groundverts.push_back(x->x);
    gtris.push_back(i - initialFaceSize);
  }
  colSys->InitGroundModel(groundverts, gtris);
  colSys->InitObjectModel(verts, otris);
#endif
}


