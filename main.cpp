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
  m.SetupBridge();
  m.SetupMouseSpring(10);

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
    m.SetMousePos(mouseX, mouseY);

    // Draw
    int pSize;
    int cSize;
    float* points = m.GetPositions2d(&pSize);
    float* colors = m.GetColors(&cSize);
    DrawDelegate::BeginFrame();
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
      printf("Frames per second: %f Springs: %d Implicit: %d\n", 1.0/secondsPerFrame, m.springs.size(), (int) implicitUpdate);
    }
  }

  glfwDestroyWindow(window);

  glfwTerminate();
  exit(EXIT_SUCCESS);
}
