#include "collision_system.h"
#include "ccdAPI.h"
#include "stdio.h"

static vec3f_list vecList;
static tri_list triList;
static std::vector<unsigned int>* vToF = NULL;
static bool initialized = false;

CollisionSystem::CollisionSystem() {}
CollisionSystem::~CollisionSystem() {
  if (initialized) {
    ccdQuitModel();
  }
}

void CollisionSystem::InitSystem(const std::vector<Eigen::Vector3d>& verts, const std::vector<int>& tris) {
  if (initialized) {
    ccdQuitModel();
  }
  initialized = true;
  vecList.clear();
  triList.clear();
  for (int i = 0; i < verts.size(); i++) {
    vecList.push_back(vec3f(verts[i][0], verts[i][1], verts[i][2]));
  }
  for (int i = 0; i < tris.size(); i += 3) {
    triList.push_back(tri3f(tris[i], tris[i + 1], tris[i + 2]));
  }
  ccdInitModel(vecList, triList);
}

void CollisionSystem::UpdateVertex(unsigned int index, const Eigen::Vector3d& vec) {
  vecList[index].x = vec[0];
  vecList[index].y = vec[1];
  vecList[index].z = vec[2];
}
void EECallback(unsigned int e1_v1, unsigned e1_v2,
				unsigned int e2_v1, unsigned int e2_v2, float t) {
	//printf("EE result: e1(%d, %d), e2(%d, %d) @ t=%f\n", e1_v1, e1_v2, e2_v1, e2_v2, t);
}
void VFCallback(unsigned int vid, unsigned int fid, float t) {
	//printf("VF result: v=%d, f=%d @ t=%f\n", vid, fid, t);
  vToF->push_back(vid);
  vToF->push_back(fid);
}

void CollisionSystem::GetCollisions(std::vector<unsigned int>& vertexToFace) {
  vToF = &vertexToFace;
  ccdUpdateVtxs(vecList);
  ccdSetEECallback(EECallback);
  ccdSetVFCallback(VFCallback);
  ccdChecking(true);
}
