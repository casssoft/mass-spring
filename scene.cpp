#include "scene.h"

#include <GLFW/glfw3.h>
#include "draw_delegate.h"

#include <stdlib.h>
#include <chrono> // milliseconds
#include <thread> // this_thread::sleep_for
#include <unistd.h> // usleep

#include <cmath>
#include <iostream>
#include "particle_system.h"

#include "Eigen/Dense"
#include "Eigen/LU"

Scene::Scene() {
  limitFps = true;
  xpos = ypos = 0;
  zpos = -20;
  xtarg = ytarg = 0;
  ztarg = 0;
  walkUp = walkDown = walkForward = walkBack = walkRight = walkLeft = false;
  drawMode = 0;
  slowMode = false;
}

void Scene::InitTime() {
  timestart = glfwGetTime();
  curTime = timestart;
  frames = 0;
  timeAhead = 0;
  prevTime = timestart;
  secondsPerFrame = 1.0/60.0;
}

void Scene::EndOfFrame() {
  frames++;
  prevTime = curTime;
  curTime = glfwGetTime();

  if (limitFps) {
    timeAhead -= curTime - prevTime;
    timeAhead += 1.0/60.0;
    if (timeAhead > 0) {
#ifndef MACOSX
      usleep((long)(timeAhead*1000000));
#else
      std::this_thread::sleep_for(std::chrono::duration<double>(timeAhead));
#endif
    }
  }
  if (curTime != timestart) {
    secondsPerFrame = (curTime - timestart)/frames;
  }
  if (secondsPerFrame * frames > .5) {
    fpsVec.push_back(1.0/secondsPerFrame);
    RecalculateFps();
  }
  if (secondsPerFrame * frames > 4) {
    printf("Frames per second: %f\n", 1.0/secondsPerFrame);
  }
}

void Scene::RecalculateFps() {
  timestart = curTime;
  frames = 0;
  timeAhead = 0;
}

double Scene::GetTimestep() {
  if (slowMode) return secondsPerFrame/20;
  return secondsPerFrame;
}

void Scene::SetLimitFps(bool enable) {
  limitFps = enable;
  RecalculateFps();
}

void Scene::ToggleLimitFps() {
  if (limitFps) limitFps = false;
  else limitFps = true;
  RecalculateFps();
}

double Scene::GetFps() {
  if (secondsPerFrame != 0)
    return 1/secondsPerFrame;
  else return 60;
}
namespace {
void RotationMatrix(const Eigen::Vector3f& position, const Eigen::Vector3f& target, const Eigen::Vector3f& up, Eigen::Matrix4f& viewMatrix) {
  Eigen::Matrix3f R;
  R.col(2) = (position-target).normalized();
  R.col(0) = up.cross(R.col(2)).normalized();
  R.col(1) = R.col(2).cross(R.col(0));
  viewMatrix.topLeftCorner<3,3>() = R.transpose();
  viewMatrix.topRightCorner<3,1>() = -R.transpose() * position;
  viewMatrix(3,3) = 1.0f;
}

void PerspectiveMatrix(float fovY, float aspect, float near, float far, Eigen::Matrix4f& projectionMatrix) {
  float theta = fovY*0.5;
  float range = far - near;
  float invtan = 1./tan(theta);

  projectionMatrix(0,0) = invtan / aspect;
  projectionMatrix(1,1) = invtan;
  projectionMatrix(2,2) = -(near + far) / range;
  projectionMatrix(2,3) = -1;
  projectionMatrix(3,2) = -2 * near * far / range;
  projectionMatrix(3,3) = 0;
}

#define PI 3.14159265
}
static Eigen::Matrix4f g_viewMatrix;
void Scene::DrawScene(ParticleSystem* m, double strainSize, bool drawPoints) {
  int pSize;
  int cSize;
  //float* points = m->GetPositions3d(&pSize);
  //float* colors = m->GetColors(&cSize, strainSize);

  float *points, *colors;
  switch(drawMode) {
    case 0:
      points = m->GetSurfaceTriangles3d(&pSize);
      colors = m->GetTriColors(&cSize, strainSize);
      break;
    case 1:
      points = m->GetPositions3d(&pSize);
      colors = m->GetColors(&cSize, strainSize, xpos, ypos, zpos);
      break;
    case 2:
      points = m->GetSurfaceTriangles3d(&pSize);
      colors = m->GetStrainSurfaceTriColors(&cSize, strainSize);
      break;
    case 3:
      points = m->GetTetCenter(&pSize);
      colors = m->GetCenterColors(&cSize, strainSize);
      break;
    case 4:
      points = m->GetAllTriangles3d(&pSize);
      colors = m->GetStrainAllTriColors(&cSize, strainSize);
      break;
  }

  Eigen::Matrix4f rotationMatrix;
  Eigen::Matrix4f projectionMatrix;
  rotationMatrix.setZero();
  projectionMatrix.setZero();
  Eigen::Vector3f pos, target, up;
  pos << xpos, ypos, zpos;
  m->GetCameraPosAndSize(&xtarg, &ytarg, &ztarg);
  target << xtarg, ytarg, ztarg;
  up << 0, -1, 0;
  RotationMatrix(pos, target, up, rotationMatrix);
  PerspectiveMatrix((65*PI)/180.0, ((float)DDWIDTH)/DDHEIGHT, .5, 100, projectionMatrix);
  g_viewMatrix = projectionMatrix * rotationMatrix;


  DrawDelegate::BeginFrame();
  DrawDelegate::SetViewMatrix(g_viewMatrix.data());
  Scene::DrawGrid(1);
  if (drawPoints) {
     DrawDelegate::SetLineSize(3);

     switch(drawMode) {
       case 0:
         DrawDelegate::DrawTriangles(points, pSize, colors, cSize);
         break;
       case 1:
         DrawDelegate::DrawLines(points, pSize, colors, cSize);
         break;
       case 2:
         DrawDelegate::DrawTriangles(points, pSize, colors, cSize);
         break;
       case 3:
         DrawDelegate::DrawPoints(points, pSize, colors, cSize);
         break;
       case 4:
         DrawDelegate::DrawTriangles(points, pSize, colors, cSize);
         break;
     }
  }
}

void Scene::Update(double timestep) {
  Eigen::Vector3d pos, targ, up, walkvector;
  pos << xpos, ypos, zpos;
  targ << xtarg, ytarg, ztarg;
  targ = targ - pos;
  targ.normalize();
  up << 0, -1, 0;
  walkvector << 0, 0, 0;
  if (walkUp) walkvector += up;
  if (walkDown) walkvector -= up;
  if (walkForward) walkvector += targ;
  if (walkBack) walkvector -= targ;
  if (walkLeft) walkvector += up.cross(targ);
  if (walkRight) walkvector += targ.cross(up);
  if (walkvector.norm() > 0)
    walkvector.normalize();
  pos += walkvector * timestep * 5;
  xpos = pos[0];
  ypos = pos[1];
  zpos = pos[2];

}
void Scene::DrawGrid(int gridSize) {
  int x_flr = -20;
  int y_flr = -20;
  int xGrids = 40;
  int yGrids = 40;
  gridpoints.resize((xGrids + yGrids) * 6);
  gridcolors.resize((xGrids + yGrids) * 6);
  for (int i = 0; i < xGrids; ++i) {
    gridpoints[i*6] = (i*gridSize + x_flr);
    gridpoints[i*6 + 1] = y_flr;
    gridpoints[i*6 + 2] = 3;
    gridpoints[i*6 + 3] = (i*gridSize + x_flr);
    gridpoints[i*6 + 4] = DDHEIGHT;
    gridpoints[i*6 + 5] = 3;
    gridcolors[i*6] = 0;
    gridcolors[i*6 + 1] = 1;
    gridcolors[i*6 + 2] = 0;
    gridcolors[i*6 + 3] = 0;
    gridcolors[i*6 + 4] = 1;
    gridcolors[i*6 + 5] = 0;
  }
  for (int i = xGrids; i < xGrids + yGrids; ++i) {
    gridpoints[i*6] = x_flr;
    gridpoints[i*6 + 1] = ((i - xGrids) *gridSize + y_flr);
    gridpoints[i*6 + 2] = 3;
    gridpoints[i*6 + 3] = DDWIDTH;
    gridpoints[i*6 + 4] = ((i - xGrids)* gridSize + y_flr);
    gridpoints[i*6 + 5] = 3;
    gridcolors[i*6] = 0;
    gridcolors[i*6 + 1] = 1;
    gridcolors[i*6 + 2] = 0;
    gridcolors[i*6 + 3] = 0;
    gridcolors[i*6 + 4] = 1;
    gridcolors[i*6 + 5] = 0;
  }
  DrawDelegate::SetLineSize(1);
  DrawDelegate::SetLineSize(1);
  DrawDelegate::DrawLines(gridpoints.data(), gridpoints.size(), gridcolors.data(), gridcolors.size());
}

int Scene::GridFloor(double x, int mpG) {
  int ret = (int)(x/mpG);
  if (x <0) ret -= 1;
  return ret * mpG;
}

void Scene::GetCameraRay(double x, double y, Eigen::Vector3d* origin, Eigen::Vector3d* ray) {
  Eigen::Matrix4f inverse;
  inverse = g_viewMatrix.inverse();
  Eigen::Vector4f preVec;
  preVec << (2 * x / DDWIDTH) - 1, 1 - (2 * y / DDHEIGHT), 2 * .5 - 1, 1;
  (*origin)[0] = xpos;
  (*origin)[1] = ypos;
  (*origin)[2] = zpos;
  Eigen::Vector4f ori = inverse * preVec;
  (*ray)[0] = ori[0] - (*origin)[0];
  (*ray)[1] = ori[1] - (*origin)[1];
  (*ray)[2] = ori[2] - (*origin)[2];
}
