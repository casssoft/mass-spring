#include "meshgen.h"
#define TETLIBRARY
#include "tetgen.h"
#include "Eigen/Sparse"
#include <vector>

void MeshGen::GenerateBar(double*& points, int& psize, std::vector<int>& tets, std::vector<int>& faces, std::vector<int>& facetotet) {
  tetgenio in, out;
  tetgenio::facet *f;
  tetgenio::polygon *p;
  int i;

  in.firstnumber = 1;

  in.numberofpoints = 8;
  in.pointlist = new REAL[in.numberofpoints * 3];
  in.pointlist[0]  = 0;  // node 1.
  in.pointlist[1]  = 0;
  in.pointlist[2]  = 0;
  in.pointlist[3]  = 2;  // node 2.
  in.pointlist[4]  = 0;
  in.pointlist[5]  = 0;
  in.pointlist[6]  = 2;  // node 3.
  in.pointlist[7]  = 2;
  in.pointlist[8]  = 0;
  in.pointlist[9]  = 0;  // node 4.
  in.pointlist[10] = 2;
  in.pointlist[11] = 0;
  // Set node 5, 6, 7, 8.
  for (i = 4; i < 8; i++) {
    in.pointlist[i * 3]     = in.pointlist[(i - 4) * 3];
    in.pointlist[i * 3 + 1] = in.pointlist[(i - 4) * 3 + 1];
    in.pointlist[i * 3 + 2] = 10;
  }

  in.numberoffacets = 6;
  in.facetlist = new tetgenio::facet[in.numberoffacets];
  in.facetmarkerlist = new int[in.numberoffacets];

  // Facet 1. The leftmost facet.
  f = &in.facetlist[0];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 1;
  p->vertexlist[1] = 2;
  p->vertexlist[2] = 3;
  p->vertexlist[3] = 4;
  
  // Facet 2. The rightmost facet.
  f = &in.facetlist[1];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 5;
  p->vertexlist[1] = 6;
  p->vertexlist[2] = 7;
  p->vertexlist[3] = 8;

  // Facet 3. The bottom facet.
  f = &in.facetlist[2];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 1;
  p->vertexlist[1] = 5;
  p->vertexlist[2] = 6;
  p->vertexlist[3] = 2;

  // Facet 4. The back facet.
  f = &in.facetlist[3];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 2;
  p->vertexlist[1] = 6;
  p->vertexlist[2] = 7;
  p->vertexlist[3] = 3;

  // Facet 5. The top facet.
  f = &in.facetlist[4];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 3;
  p->vertexlist[1] = 7;
  p->vertexlist[2] = 8;
  p->vertexlist[3] = 4;

  // Facet 6. The front facet.
  f = &in.facetlist[5];
  f->numberofpolygons = 1;
  f->polygonlist = new tetgenio::polygon[f->numberofpolygons];
  f->numberofholes = 0;
  f->holelist = NULL;
  p = &f->polygonlist[0];
  p->numberofvertices = 4;
  p->vertexlist = new int[p->numberofvertices];
  p->vertexlist[0] = 4;
  p->vertexlist[1] = 8;
  p->vertexlist[2] = 5;
  p->vertexlist[3] = 1;

  // Set 'in.facetmarkerlist'

  in.facetmarkerlist[0] = -1;
  in.facetmarkerlist[1] = -2;
  in.facetmarkerlist[2] = 0;
  in.facetmarkerlist[3] = 0;
  in.facetmarkerlist[4] = 0;
  in.facetmarkerlist[5] = 0;

  // Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
  //   do quality mesh generation (q) with a specified quality bound
  //   (1.414), and apply a maximum volume constraint (a0.1).

  tetrahedralize("pnnq1.414a.5", &in, &out);

  points = new double[out.numberofpoints*3];
  for (int i = 0; i < out.numberofpoints*3;++i) {
    if (i%3 == 2) {
      points[i] = -1* out.pointlist[i];
    } else {
      points[i] = out.pointlist[i];
    }
  }
  psize = out.numberofpoints;

  for (int i = 0; i <out.numberoftetrahedra*4; ++i) {
    tets.push_back(out.tetrahedronlist[i] - 1);
  }
  for (int i = 0; i < out.numberoftrifaces*3; ++i) {
    faces.push_back(out.trifacelist[i] - 1);
  }
  for (int i = 0; i < out.numberoftrifaces; ++i) {
    if (out.adjtetlist[i*2] - 1< 0 || out.adjtetlist[i*2] - 1 >= out.numberoftetrahedra) {
      if (out.adjtetlist[i*2+1] - 1< 0 || out.adjtetlist[i*2+1] - 1 >= out.numberoftetrahedra) {
        printf("No adj tet for this face %d\n", i);
        facetotet.push_back(0);
      } else {
        facetotet.push_back(out.adjtetlist[i*2+1] - 1);
      }
    } else {
      facetotet.push_back(out.adjtetlist[i*2] - 1);
    }
  }
}

void MeshGen::GenerateMesh(double*& points, int& psize, std::vector<int>& tets, std::vector<int>& faces, std::vector<int>& facetotet, char* filename) {
  tetgenio out, in;
  int i;

  // All indices start from 0.
  in.firstnumber = 0;
  in.load_ply(filename);

  // Tetrahedralize the PLC. Switches are chosen to read a PLC (p),
  //   do quality mesh generation (q) with a specified quality bound
  //   (1.414), and apply a maximum volume constraint (a0.1).

  tetrahedralize("pnnqa.05", &in, &out);

  points = new double[out.numberofpoints*3];
  for (int i = 0; i < out.numberofpoints*3;++i) {
    if (i%3 == 2) {
      points[i] = -1* out.pointlist[i];
    } else {
      points[i] = out.pointlist[i];
    }
  }
  psize = out.numberofpoints;

  for (int i = 0; i <out.numberoftetrahedra*4; ++i) {
    tets.push_back(out.tetrahedronlist[i]);
  }
  for (int i = 0; i < out.numberoftrifaces*3; ++i) {
    faces.push_back(out.trifacelist[i]);
  }
  for (int i = 0; i < out.numberoftrifaces; ++i) {
    if (out.adjtetlist[i*2] < 0 || out.adjtetlist[i*2] >= out.numberoftetrahedra) {
      if (out.adjtetlist[i*2+1] < 0 || out.adjtetlist[i*2+1] >= out.numberoftetrahedra) {
        facetotet.push_back(0);
      } else {
        facetotet.push_back(out.adjtetlist[i*2+1]);
      }
    } else {
      facetotet.push_back(out.adjtetlist[i*2]);
    }
  }
}
