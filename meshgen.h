#ifndef MESHGEN_H__
#define MESHGEN_H__

#include <vector>
namespace MeshGen {
  void GenerateBar(double*& points, int& psize, std::vector<int>& edges);
  void GenerateMesh(double*& points, int& psize, std::vector<int>& edges, std::vector<int>& faces, char*filename);
};

#endif
