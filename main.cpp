#include <GLFW/glfw3.h>

#include "draw_delegate.h"
#include "particle_system.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "chrono" // milliseconds
#include "thread" // this_thread::sleep_for
#include <unistd.h> // usleep

namespace {
void error_callback(int error, const char* description) {
  fprintf(stderr, "%s\n", description);
}
int changeSetup = 0;
bool limitFps = true;
bool recalculateFps = false;
bool implicitUpdate = false;
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (action == GLFW_PRESS) {
    switch(key) {
      case GLFW_KEY_ESCAPE:
        glfwSetWindowShouldClose(window, GL_TRUE);
        break;

      case 'Q':
        changeSetup = 1;
        break;

      case 'W':
        changeSetup = 2;
        break;

      case 'A':
        if (limitFps) limitFps = false;
        else limitFps = true;
        recalculateFps = true;
        break;

      case 'S':
        if (implicitUpdate) implicitUpdate = false;
        else implicitUpdate = true;
        break;
    }
  }
}

int GridFloor(double x, int mpG) {
  int ret = (int)(x/mpG);
  if (x <0) ret -= 1;
  return ret * mpG;
}

}; // namespace

int main(int argc, char **argv) {
  GLFWwindow* window;
  glfwSetErrorCallback(error_callback);
  if (!glfwInit())
      exit(EXIT_FAILURE);
  window = glfwCreateWindow(DDWIDTH, DDHEIGHT, "Springs", NULL, NULL);

  if (!window)
  {
      glfwTerminate();
      exit(EXIT_FAILURE);
  }

  glfwMakeContextCurrent(window);
  glfwSetKeyCallback(window, key_callback);
  glfwSwapInterval(0);
  DrawDelegate::SetupOpenGL();

  // Particle system setup
  ParticleSystem m;
  if (argc == 4) {
    m.SetSpringProperties(atof(argv[1]), atof(argv[2]));
    implicitUpdate = atoi(argv[3]);
  }

  //m.SetupSingleSpring();
  m.SetupBridge2();
  m.SetupMouseSpring(5);
  //m.SetupTriangle();

  // Frames per second set up
  double timestart = glfwGetTime();
  double curTime = timestart;
  int frames = 0;
  double timeAhead = 0;
  double prevTime = timestart;

  double secondsPerFrame = 1.0/60.0;

  // Cursor spring set up
  double mouseX, mouseY;
  glfwGetCursorPos(window, &mouseX, &mouseY);
  m.SetMousePos(mouseX, mouseY);

  // Camera setup
  double x, y, zoom;

  x = 0; y = 0; zoom = 1;
  std::vector<float> gridpoints;
  std::vector<float> gridcolors;
  while (!glfwWindowShouldClose(window)) {

    // Handle changing setup
    if (changeSetup == 1) {
      m.Reset();
      m.SetupTriforce();
      m.SetupMouseSpring(0);
      changeSetup = 0;
    } else if (changeSetup == 2) {
      m.Reset();
      m.SetupTriangle();
      m.SetupMouseSpring(0);
      changeSetup = 0;
    }

    // Handle spring enable from mouse button
    if (glfwGetMouseButton(window, 0) == GLFW_PRESS) {
      m.SetMouseSpring(false);
    } else {
      m.SetMouseSpring(true);
    }

    // Update m
    m.Update(secondsPerFrame, implicitUpdate);

    // Set mouse spring pos
    glfwGetCursorPos(window, &mouseX, &mouseY);
    m.SetMousePos(mouseX/zoom + x, mouseY/zoom + y);

    // Draw
    int pSize;
    int cSize;
    m.GetCameraPosAndSize(&x, &y, &zoom);
    float* points = m.GetPositions2d(&pSize, x, y, zoom);
    float* colors = m.GetColors(&cSize);

    // Make grid
    int mpG = 5;

    int x_flr = GridFloor(x, mpG);
    int y_flr = GridFloor(y, mpG);
    int xGrids = GridFloor(DDWIDTH/zoom + x, mpG) - GridFloor(x, mpG) + 1;
    int yGrids = GridFloor(DDHEIGHT/zoom + y, mpG) - GridFloor(y, mpG) + 1;
    gridpoints.resize((xGrids + yGrids) * 4);
    gridcolors.resize((xGrids + yGrids) * 6);
    for (int i = 0; i < xGrids; ++i) {
      gridpoints[i*4] = (i*mpG + x_flr - x) * zoom;
      gridpoints[i*4 + 1] = 0;
      gridpoints[i*4 + 2] = (i*mpG + x_flr - x) * zoom;
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
      gridpoints[i*4 + 1] = ((i - xGrids) *mpG + y_flr - y) * zoom;
      gridpoints[i*4 + 2] = DDWIDTH;
      gridpoints[i*4 + 3] = ((i - xGrids)* mpG + y_flr - y) * zoom;
      gridcolors[i*6] = 0;
      gridcolors[i*6 + 1] = 1;
      gridcolors[i*6 + 2] = 0;
      gridcolors[i*6 + 3] = 0;
      gridcolors[i*6 + 4] = 1;
      gridcolors[i*6 + 5] = 0;
    }
    DrawDelegate::BeginFrame();
    DrawDelegate::SetLineSize(1);
    DrawDelegate::DrawLines(gridpoints.data(), (xGrids + yGrids) * 4, gridcolors.data(), (xGrids + yGrids) * 6);
    DrawDelegate::SetLineSize(4);
    DrawDelegate::DrawLines(points, pSize, colors, cSize);


    glfwSwapBuffers(window);
    glfwPollEvents();


    // Handle fps
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
    if (secondsPerFrame * frames > 4 || recalculateFps) {
      recalculateFps = false;
      timestart = curTime;
      frames = 0;
      printf("Frames per second: %f Springs: %d Implicit: %d\n x: %f y:%f zoom:%f\n", 1.0/secondsPerFrame, m.springs.size(), (int) implicitUpdate, x, y, zoom);
    }
  }

  glfwDestroyWindow(window);

  glfwTerminate();
  exit(EXIT_SUCCESS);
}
