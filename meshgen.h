#ifndef MESHGEN_H__
#define MESHGEN_H__

#include <vector>
namespace MeshGen {
  void GenerateBar(double*& points, int& psize, std::vector<int>& edges, std::vector<int>& faces, std::vector<int>& facetotet);
  void GenerateMesh(double*& points, int& psize, std::vector<int>& edges, std::vector<int>& faces, std::vector<int>& facetotet, const char*filename);
};

#endif
