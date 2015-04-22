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
int bridgeL = 10;
float xpos = 0, ypos = 0, zpos = -20;
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
        zpos += 1;
        break;
      case 'K':
        zpos -= 1;
        break;
      case 'J':
        xpos -= 1;
        break;
      case 'L':
        xpos += 1;
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
  m.SetupBridge(bridgeL);
  //m.SetupTriangle();
  //m.SetupMouseSpring(5);

  Scene scene;
  scene.InitTime();
  scene_p = &scene;


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
        m.SetupBridge(bridgeL);
        break;
      case 3:
        m.SetupBendingBar();
        break;
      case 4:
        m.SetupTriforce();
        break;
    }
    changeSetup = 0;

    // Update m
    m.Update(scene.GetTimestep(), implicitUpdate);

    // Draw
    scene.DrawScene(&m, strainSize, xpos, ypos, zpos);

    glfwSwapBuffers(window);
    glfwPollEvents();

    scene.EndOfFrame();
  }

  glfwDestroyWindow(window);

  glfwTerminate();
  exit(EXIT_SUCCESS);
}
