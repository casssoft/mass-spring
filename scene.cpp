#include "scene.h"

#include <GLFW/glfw3.h>
#include "draw_delegate.h"

#include <stdlib.h>
#include <chrono> // milliseconds
#include <thread> // this_thread::sleep_for
#include <unistd.h> // usleep

#include <cmath>

#include "particle_system.h"

#include "Eigen/Dense"

Scene::Scene() {
  limitFps = true;
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
  if (secondsPerFrame * frames > 4) {
    printf("Frames per second: %f\n", 1.0/secondsPerFrame);
    RecalculateFps();
  }
}

void Scene::RecalculateFps() {
  timestart = curTime;
  frames = 0;
  timeAhead = 0;
}

double Scene::GetTimestep() {
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
void Scene::DrawScene(ParticleSystem* m, int strainSize, float xpos, float ypos, float zpos, bool drawPoints) {
  int pSize;
  int cSize;
  //float* points = m->GetPositions3d(&pSize);
  //float* colors = m->GetColors(&cSize, strainSize);

  float* points = m->GetTriangles3d(&pSize);
  float* colors = m->GetTriColors(&cSize, strainSize);
  Eigen::Matrix4f rotationMatrix;
  Eigen::Matrix4f projectionMatrix;
  rotationMatrix.setZero();
  projectionMatrix.setZero();
  Eigen::Vector3f pos, target, up;
  pos << xpos, ypos, zpos;
  target << points[0], points[1], points[2];
  up << 0, -1, 0;
  RotationMatrix(pos, target, up, rotationMatrix);
  PerspectiveMatrix((65*PI)/180.0, ((float)DDWIDTH)/DDHEIGHT, .5, 100, projectionMatrix);
  Eigen::Matrix4f viewMatrix = projectionMatrix * rotationMatrix;


  DrawDelegate::BeginFrame();
  DrawDelegate::SetViewMatrix(viewMatrix.data());
  Scene::DrawGrid(1);
  if (drawPoints) {
     DrawDelegate::SetLineSize(3);
     DrawDelegate::DrawTriangles(points, pSize, colors, cSize);
  }
  if (frames% 100 == 0) {
     printf("pSize %d\n", pSize);
  }
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
