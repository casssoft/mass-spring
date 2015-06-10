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
bool drawSimulation = true;
bool solveWithguess = true;
bool corotational = false;
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
        if (drawSimulation) drawSimulation = false;
        else drawSimulation = true;
        printf("Implicit: %d\n", (int) drawSimulation);
        break;
      case 'D':
        if (solveWithguess) solveWithguess = false;
        else solveWithguess = true;
        printf("Solvewith guess: %d\n", (int) solveWithguess);
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
      case 'U':
        scene_p->walkUp = true;
        break;
      case 'O':
        scene_p->walkDown = true;
        break;
      case 'Z':
        scene_p->drawMode = (scene_p->drawMode + 1)%5;
        break;
      case 'X':
        scene_p->slowMode = scene_p->slowMode ? false : true;
        break;
      case 'C':
        corotational = corotational ? false : true;
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
      case 'U':
        scene_p->walkUp = false;
        break;
      case 'O':
        scene_p->walkDown = false;
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

  char*meshFilename = "Armadillo_simple2.ply";
  // Particle system setup
  ParticleSystem m;
  double strainSize = 10;
  if (argc >= 3) {
    m.SetSpringProperties(atof(argv[1]), atof(argv[2]));
    if (argc >= 4) {
      strainSize = atof(argv[3]);
    }
    if (argc >= 5) {
      meshFilename = argv[4];
    }
  }

  m.SetupArmadillo();
  //m.SetupSingleSpring();
  //m.SetupBendingBar();
  //m.SetupTriangle();
  //m.SetupMouseSpring(5);

  Scene scene;
  scene.InitTime();
  scene_p = &scene;

  double curTime;
  double simulatetime = 0;
  double drawtime = 0;
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
        m.SetupArmadillo();
        break;
      case 3:
        m.SetupBendingBar();
        break;
      case 4:
        m.SetupMeshFile(meshFilename);
        break;
    }
    changeSetup = 0;
    curTime = glfwGetTime();
    // Update m
    double timestep = scene.GetTimestep();

    m.Update(timestep, solveWithguess, corotational);

    double tempTime = glfwGetTime();
    simulatetime += tempTime - curTime;
    curTime = tempTime;


    // Update scene pos
    scene.Update(timestep);

    // Draw
    scene.DrawScene(&m, strainSize, drawSimulation);

    tempTime = glfwGetTime();
    drawtime += tempTime - curTime;
    curTime = tempTime;

    glfwSwapBuffers(window);

    glfwPollEvents();
    frames++;
    scene.EndOfFrame();
  }

  printf("Simulate time %f\n", simulatetime/frames);
  printf("Draw time %f\n", drawtime/frames);
  printf("Frames %d\n", frames);
  double triplet, fromtriplet, solve;
  m.GetProfileInfo(triplet, fromtriplet, solve);
  printf("Triplet %f, from %f, solve %f \n", triplet/frames, fromtriplet/frames, solve/frames);
  printf("Total %f\n", triplet + fromtriplet + solve);
  glfwDestroyWindow(window);

  glfwTerminate();
  exit(EXIT_SUCCESS);
}
