#include <GLFW/glfw3.h>

#include "draw_delegate.h"
#include "particle_system.h"
#include "scene.h"

#include <imgui.h>
#include "imgui_impl.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <dirent.h>

#include <vector>
#include <string>
namespace {
std::vector<std::string> modelList;
bool mloaded = false;
void LoadModelList() {
  modelList.clear();
  char buffer[1000];
  DIR *dir;
  struct dirent *ent;
  if ((dir = opendir(".")) != NULL) {
    //print all the files and directories within directory
    while ((ent = readdir(dir)) != NULL) {
      if (!strcmp(ent->d_name, ".") || !strcmp(ent->d_name, "..") ) continue;
      if (ent->d_type != DT_REG) continue;
        int len = strlen(ent->d_name);
        if (strncmp(ent->d_name + len - 4, ".ply", 5) != 0) {
          //printf("Found %s not ply file\n", ent->d_name);
          continue;
        }
        modelList.push_back(ent->d_name);
    }
    closedir(dir);
    mloaded = true;
    return;
  }
}

void error_callback(int error, const char* description) {
  fprintf(stderr, "%s\n", description);
}
bool drawSimulation = true;
bool solveWithguess = true;
bool corotational = true;
int typeOfGround = 5;
Scene* scene_p;
void key_callback(GLFWwindow* window, int key, int scancode, int action, int mods)
{
  if (action == GLFW_PRESS) {
    switch(key) {
      case GLFW_KEY_ESCAPE:
        glfwSetWindowShouldClose(window, GL_TRUE);
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
      case 'G':
        typeOfGround++;
        printf("Type of ground : %i\n", typeOfGround);
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

int main(int argc, const char **argv) {
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
  ImGui_ImplGlfw_Init(window, true);


  const char*meshFilename = "Armadillo_simple2.ply";
  // Particle system setup
  ParticleSystem m;
  double strainSize = 2;
  if (argc >= 4) {
    m.SetSpringProperties(atof(argv[1]), .4, atof(argv[2]), 9.8f, atof(argv[3]));
    if (argc >= 5) {
      strainSize = atof(argv[4]);
    }
    if (argc >= 6) {
      meshFilename = argv[5];
    }
  }

  //m.SetupArmadillo();
  //m.SetupSingleSpring();
  m.SetupBendingBar();
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

    glfwPollEvents();
    ImGui_ImplGlfw_NewFrame();

    curTime = glfwGetTime();
    // Update m
    double timestep = scene.GetTimestep();

    m.Update(timestep, solveWithguess, corotational, typeOfGround);

    double tempTime = glfwGetTime();
    simulatetime += tempTime - curTime;
    curTime = tempTime;


    // Update scene pos
    scene.Update(timestep);

    // Draw
    scene.DrawScene(&m, strainSize, drawSimulation);

    {
      ImGui::Text("Change configuration");
      static int selected_config = 2;
      const char* names[] = { "Single Spring", "Armadillo", "Bending Bar", "Load Mesh File" };
      int nameLength = 4;
      const char* groundTypes[] = {
        "No Ground",
        "Only Penalty",
        "Snap to floor and penalty",
        "Snap to prev intersection and penalty",
        "Snap to prev intersection and penalty plus friction penalty",
        "Snap to floor and infinite friction",
        "Snap to floor and implicit penalty"
      };
      int groundTypeLength = 7;

      if (ImGui::Button("Select Type.."))
          ImGui::OpenPopup("select");
      ImGui::SameLine();
      ImGui::Text(selected_config == -1 ? "<None>" : names[selected_config]);
      if (ImGui::BeginPopup("select"))
      {
          for (int i = 0; i < 4; i++)
              if (ImGui::Selectable(names[i]))
                  selected_config = i;
          ImGui::EndPopup();
      }
      if (selected_config == 3) {
        if (ImGui::Button("Select File.."))
            ImGui::OpenPopup("select_file");
        ImGui::SameLine();
        ImGui::Text((meshFilename != NULL && meshFilename[0] != 0) ? meshFilename : "<None>");

        if (ImGui::BeginPopup("select_file")) {
          if (!mloaded) LoadModelList();
          for (int i = 0; i < modelList.size(); i++) {
              if (ImGui::Selectable(modelList[i].c_str())) {
                  meshFilename = modelList[i].c_str();
              }
          }
          ImGui::Separator();
          if (ImGui::Selectable("Reload List")) {
            LoadModelList();
            meshFilename = NULL;
          }
          ImGui::EndPopup();
        }
      }

      ImGui::Separator();

      if (ImGui::Button("Select Ground Type.."))
          ImGui::OpenPopup("select_ground");
      ImGui::SameLine();
      ImGui::Text(typeOfGround < 0 && typeOfGround >= groundTypeLength ? "Undefined ground" : groundTypes[typeOfGround]);
      if (ImGui::BeginPopup("select_ground"))
      {
          for (int i = 0; i < groundTypeLength; i++)
              if (ImGui::Selectable(groundTypes[i]))
                  typeOfGround = i;
          ImGui::EndPopup();
      }

      ImGui::Separator();

      const char* drawModeTypes[] = {
        "Normal",
        "Lines between tets",
        "Surface triangles and strain (broken for load mesh)",
        "Tet center with strain",
        "All tets and strain"
      };
      int drawModeLength = 5;
      if (ImGui::Button("Select Draw Mode.."))
          ImGui::OpenPopup("select_drawMode");
      ImGui::SameLine();
      ImGui::Text(scene_p->drawMode < 0 && scene_p->drawMode >= drawModeLength ? "Undefined draw mode" : drawModeTypes[scene_p->drawMode]);
      if (ImGui::BeginPopup("select_drawMode"))
      {
          for (int i = 0; i < drawModeLength; i++)
              if (ImGui::Selectable(drawModeTypes[i]))
                  scene_p->drawMode = i;
          ImGui::EndPopup();
      }

      ImGui::Separator();

      static float stiffness = 1000.0f;
      static float volumeConservation = 0.4f;
      static float damping = 10.0f;
      static float gravity = 9.6f;
      static float strainDisplaySize = strainSize;
      static float groundStiffness = 1000.0f;
      ImGui::Text("Stiffness (Young's modulus)");
      ImGui::SliderFloat("##stiffness", &stiffness, 0.0f, 10000.0f, "%.3f", 2.0);
      ImGui::Text("Volume Conservation (Poisson's ratio)");
      ImGui::SliderFloat("##vol", &volumeConservation, 0.0f, .49999f);
      ImGui::Text("Damping (Not supported with FEM yet)");
      ImGui::SliderFloat("##damp", &damping, 0.0f, 1000.0f, "%.3f", 2.0);
      ImGui::Text("Gravity m/s^2");
      ImGui::SliderFloat("##grav", &gravity, 0.0f, 50.0f);
      ImGui::Text("Strain display ratio");
      ImGui::SliderFloat("##strainSize", &strainDisplaySize, 0.0f, 100.0f, "%.3f", 4.0);
      ImGui::Text("Ground stiffness (size of penalty forces)");
      ImGui::SliderFloat("##groundstiffness", &groundStiffness, 0.0f, 10000.0f);

      if (ImGui::Button("Apply Changes")) {
        m.SetSpringProperties(stiffness, volumeConservation, damping, gravity, groundStiffness);
        strainSize = strainDisplaySize;
        switch (selected_config) {
          case 0:
            m.SetupSingleSpring();
            break;
          case 1:
            m.SetupArmadillo();
            break;
          case 2:
            m.SetupBendingBar();
            break;
          case 3:
            if (meshFilename != NULL && meshFilename[0] != 0)
              m.SetupMeshFile(meshFilename);
            break;
        }
      }
    }
    if (scene_p->fpsVec.size() != 0) {
      int windowSize = 20;
      int offset = 0;
      if (scene_p->fpsVec.size() > 40) {
        offset = scene_p->fpsVec.size() - 40;
      }
      ImGui::PlotLines("Frames Per Second##Graph", scene_p->fpsVec.data() + offset, scene_p->fpsVec.size() - offset, 0, NULL, 0.0, 200.0f, ImVec2(0,150));
      char buffer[1000];
      snprintf(buffer, 1000, "Current FPS: %.3f", scene_p->fpsVec[scene_p->fpsVec.size() - 1]);
      ImGui::Text(buffer);
      ImGui::Checkbox("Limit to 60 FPS?", &(scene_p->limitFps));
    }

    ImGui::Render();

    tempTime = glfwGetTime();
    drawtime += tempTime - curTime;
    curTime = tempTime;

    glfwSwapBuffers(window);

    frames++;
    scene.EndOfFrame();
  }

  printf("Simulate time %f\n", simulatetime/frames);
  printf("Draw time %f\n", drawtime/frames);
  printf("Frames %d\n", frames);
  double triplet, fromtriplet, solve, setup;
  m.GetProfileInfo(triplet, fromtriplet, solve, setup);
  printf("Triplet %f, from %f, solve %f, setup %f\n", triplet/frames, fromtriplet/frames, solve/frames, setup/frames);
  printf("Total %f\n", triplet + fromtriplet + solve + setup);

  ImGui_ImplGlfw_Shutdown();

  glfwDestroyWindow(window);

  glfwTerminate();
  exit(EXIT_SUCCESS);
}
