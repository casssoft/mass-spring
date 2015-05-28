#include <GLFW/glfw3.h>

#include "draw_delegate.h"
#include "particle_system.h"
#include "scene.h"

#include "stdlib.h"
#include "stdio.h"
#include "string.h"


namespace {
void error_callback(int error, const char* description) {
  fprintf(stderr, "%s\n", description);
}
int changeSetup = 0;
bool implicitUpdate = false;
bool solveWithguess = true;
int bridgeL = 10;
Scene* scene_p;
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (action == GLFW_PRESS) {
    switch(key) {
      case GLFW_KEY_ESCAPE:
        glfwSetWindowShouldClose(window, GL_TRUE);
        break;

      case '1':
        changeSetup = 1;
        break;

      case '2':
        changeSetup = 2;
        break;

      case '3':
        changeSetup = 3;
        break;

      case '4':
        changeSetup = 4;
        break;

      case 'A':
        scene_p->ToggleLimitFps();
        break;

      case 'S':
        if (implicitUpdate) implicitUpdate = false;
        else implicitUpdate = true;
        printf("Implicit: %d\n", (int) implicitUpdate);
        break;

      case 'Q':
        bridgeL++;
        break;

      case 'W':
        bridgeL--;
        break;
      case 'I':
        scene_p->walkForward = true;
        break;
      case 'K':
        scene_p->walkBack = true;
        break;
      case 'J':
        scene_p->walkLeft = true;
        break;
      case 'L':
        scene_p->walkRight = true;
        break;
      case 'D':
        if (solveWithguess) solveWithguess = false;
        else solveWithguess = true;
        printf("Solvewith guess: %d\n", (int) solveWithguess);
      case 'Z':
        scene_p->displaySurface = scene_p->displaySurface ? false : true;
        break;
      case 'X':
        scene_p->slowMode = scene_p->slowMode ? false : true;
        break;
    }
  } else if (action == GLFW_RELEASE) {
    switch(key) {
      case 'I':
        scene_p->walkForward = false;
        break;
      case 'K':
        scene_p->walkBack = false;
        break;
      case 'J':
        scene_p->walkLeft = false;
        break;
      case 'L':
        scene_p->walkRight = false;
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
  glfwSwapInterval(0);
  DrawDelegate::SetupOpenGL();

  // Particle system setup
  ParticleSystem m;
  int strainSize = 10;
  if (argc >= 4) {
    m.SetSpringProperties(atof(argv[1]), atof(argv[2]));
    implicitUpdate = atoi(argv[3]);
    if (argc == 5) {
      strainSize = atoi(argv[4]);
    }
  }

  //m.SetupSingleSpring();
  m.SetupBendingBar();
  //m.SetupTriangle();
  //m.SetupMouseSpring(5);

  Scene scene;
  scene.InitTime();
  scene_p = &scene;

  double curTime;
  double simulatetime = 0;
  int frames = 0;
  glfwSetKeyCallback(window, key_callback);
  while (!glfwWindowShouldClose(window)) {

    // Handle changing setup
    switch (changeSetup) {
      case 0:
        break;
      case 1:
        m.SetupSingleSpring();
        break;
      case 2:
        //m.SetupBridge(bridgeL);
        break;
      case 3:
        m.SetupBendingBar();
        break;
      case 4:
        //m.SetupTriforce();
        break;
    }
    changeSetup = 0;
    curTime = glfwGetTime();
    // Update m
    double timestep = scene.GetTimestep();

    m.Update(timestep, implicitUpdate, solveWithguess);

    double tempTime = glfwGetTime();
    simulatetime += tempTime - curTime;
    curTime = tempTime;

    // Update scene pos
    scene.Update(timestep);

    // Draw
    scene.DrawScene(&m, strainSize, implicitUpdate);

    glfwSwapBuffers(window);

    glfwPollEvents();
    frames++;
    scene.EndOfFrame();
  }

  printf("Simulate time %f\n", simulatetime/frames);
  printf("Frames %d\n", frames);
  double triplet, fromtriplet, solve;
  m.GetProfileInfo(triplet, fromtriplet, solve);
  printf("Triplet %f, from %f, solve %f \n", triplet/frames, fromtriplet/frames, solve/frames);
  printf("Total %f\n", triplet + fromtriplet + solve);
  glfwDestroyWindow(window);

  glfwTerminate();
  exit(EXIT_SUCCESS);
}
