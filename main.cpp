#include <GLFW/glfw3.h>

#include "draw_delegate.h"
#include "particle_system.h"
#include "stdlib.h"
#include "stdio.h"
#include "string.h"

#include "chrono" // milliseconds
#include "thread" // this_thread::sleep_for

namespace {
void error_callback(int error, const char* description) {
  fprintf(stderr, "%s\n", description);
}
}

int main(void)
{
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

  glfwSwapInterval(0);
  DrawDelegate::SetupOpenGL();

  // Particle system setup
  ParticleSystem m;
  m.SetupTriangle();
  m.SetupMouseSpring(0);

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

    if (glfwGetMouseButton(window, 0) == GLFW_PRESS) {
      m.SetMouseSpring(false);
    } else {
      m.SetMouseSpring(true);
    }
    m.Update(secondsPerFrame);

    glfwGetCursorPos(window, &mouseX, &mouseY);
    m.SetMousePos(mouseX, mouseY);

    int pSize;
    float* points = m.GetPositions2d(&pSize);

    DrawDelegate::BeginFrame();
    DrawDelegate::DrawLines(points, pSize);
    glfwSwapBuffers(window);
    glfwPollEvents();


    frames++;
    prevTime = curTime;
    curTime = glfwGetTime();

    timeAhead -= curTime - prevTime;
    timeAhead += 1.0/60.0;
    if (timeAhead > 0) {
      //std::this_thread::sleep_for(std::chrono::duration<double>(timeAhead));
    }
    if (curTime != timestart) {
      secondsPerFrame = (curTime - timestart)/frames;
    }
    if (secondsPerFrame * frames > 5) {
      timestart = curTime;
      frames = 0;
      printf("Frames per second: %f Springs: %d\n", 1.0/secondsPerFrame, m.springs.size());
    }
  }

  glfwDestroyWindow(window);

  glfwTerminate();
  exit(EXIT_SUCCESS);
}
