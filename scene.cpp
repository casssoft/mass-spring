#include "scene.h"

#include <GLFW/glfw3.h>
#include "draw_delegate.h"

#include <stdlib.h>
#include <chrono> // milliseconds
#include <thread> // this_thread::sleep_for
#include <unistd.h> // usleep


Scene::Scene() {
  cam_x = cam_y = 0;
  zoom = 0;
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
    printf("Frames per second: %f x: %f y:%f zoom:%f\n", 1.0/secondsPerFrame, cam_x, cam_y, zoom);
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

void Scene::DrawGrid(int gridSize) {
  int x_flr = GridFloor(cam_x, gridSize);
  int y_flr = GridFloor(cam_y, gridSize);
  int xGrids = GridFloor(DDWIDTH/zoom + cam_x, gridSize) - GridFloor(cam_x, gridSize) + 1;
  int yGrids = GridFloor(DDHEIGHT/zoom + cam_y, gridSize) - GridFloor(cam_y, gridSize) + 1;
  gridpoints.resize((xGrids + yGrids) * 4);
  gridcolors.resize((xGrids + yGrids) * 6);
  for (int i = 0; i < xGrids; ++i) {
    gridpoints[i*4] = (i*gridSize + x_flr - cam_x) * zoom;
    gridpoints[i*4 + 1] = 0;
    gridpoints[i*4 + 2] = (i*gridSize + x_flr - cam_x) * zoom;
    gridpoints[i*4 + 3] = DDHEIGHT;
    gridcolors[i*6] = 0;
    gridcolors[i*6 + 1] = 1;
    gridcolors[i*6 + 2] = 0;
    gridcolors[i*6 + 3] = 0;
    gridcolors[i*6 + 4] = 1;
    gridcolors[i*6 + 5] = 0;
  }
  for (int i = xGrids; i < xGrids + yGrids; ++i) {
    gridpoints[i*4] = 0;
    gridpoints[i*4 + 1] = ((i - xGrids) *gridSize + y_flr - cam_y) * zoom;
    gridpoints[i*4 + 2] = DDWIDTH;
    gridpoints[i*4 + 3] = ((i - xGrids)* gridSize + y_flr - cam_y) * zoom;
    gridcolors[i*6] = 0;
    gridcolors[i*6 + 1] = 1;
    gridcolors[i*6 + 2] = 0;
    gridcolors[i*6 + 3] = 0;
    gridcolors[i*6 + 4] = 1;
    gridcolors[i*6 + 5] = 0;
  }
  DrawDelegate::SetLineSize(1);
  DrawDelegate::DrawLines(gridpoints.data(), gridpoints.size(), gridcolors.data(), gridcolors.size());
}

int Scene::GridFloor(double x, int mpG) {
  int ret = (int)(x/mpG);
  if (x <0) ret -= 1;
  return ret * mpG;
}
