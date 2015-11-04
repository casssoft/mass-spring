#ifndef COLLISION_SYSTEM_PQP_H__
#define COLLISION_SYSTEM_PQP_H__
#include "../Eigen/Core"
#include <vector>

class CollisionSystemPQP {
 public:
  CollisionSystemPQP();
  ~CollisionSystemPQP();
  void InitGroundModel(const std::vector<Eigen::Vector3d>& verts, const std::vector<int>& tris);
  void InitObjectModel(const std::vector<Eigen::Vector3d>& verts, const std::vector<int>& tris);

  void GetCollisions(std::vector<unsigned int>& objectVertexToFace, std::vector<unsigned int>& edgeToEdge, std::vector<double>& edgeU, std::vector<Eigen::Vector3d>& moveEdge);
};
#endif
