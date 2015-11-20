#include "collision_system_pqp.h"
#include <Eigen/Dense>
#include <math.h>
#define FABS(x) (double(fabs(x)))        /* implement as is fastest on your machine */
#include "stdio.h"
#include "PQP.h"
static bool initialized = false;

static bool vflip, uflip;
static double vinter[2];
static double uinter[2];
static int tri_ve1[2], tri_ve2[2];
static int tri_ue1[2], tri_ue2[2];
CollisionSystemPQP::CollisionSystemPQP() {}
CollisionSystemPQP::~CollisionSystemPQP() {
  if (initialized) {
  }
}

// first model is ground the rest can move
//std::vector<PQP_Model> models;
PQP_Model ground;
PQP_Model object;
static int NoDivTriTriIsect(double V0[3],double V1[3],double V2[3],
                     double U0[3],double U1[3],double U2[3]);
static int 
TriContact(PQP_REAL *P1, PQP_REAL *P2, PQP_REAL *P3,
           PQP_REAL *Q1, PQP_REAL *Q2, PQP_REAL *Q3);

static std::vector<Eigen::Vector3d> overts;
static std::vector<Eigen::Vector3d> gverts;
static std::vector<int> otris;
static std::vector<int> gtris;
static void AddToModel(PQP_Model& m,const std::vector<Eigen::Vector3d>& verts, const std::vector<int>& tris) {
  PQP_REAL p1[3], p2[3], p3[3];
  for (int i = 0; i < tris.size(); i += 3) {
    p1[0] = verts[tris[i]][0];
    p1[1] = verts[tris[i]][1];
    p1[2] = verts[tris[i]][2];
    
    p2[0] = verts[tris[i + 1]][0];
    p2[1] = verts[tris[i + 1]][1];
    p2[2] = verts[tris[i + 1]][2];
    
    p3[0] = verts[tris[i + 2]][0];
    p3[1] = verts[tris[i + 2]][1];
    p3[2] = verts[tris[i + 2]][2];
    //fprintf(stderr, "i: %i i/3 %i\n", i, i/3);
    m.AddTri(p1, p2, p3, i/3);
  }
}

void CollisionSystemPQP::InitGroundModel(const std::vector<Eigen::Vector3d>& verts, const std::vector<int>& tris) {
  gverts = verts;
  gtris = tris;
  ground.BeginModel();
  AddToModel(ground, verts, tris);
  ground.EndModel();
}

void CollisionSystemPQP::InitObjectModel(const std::vector<Eigen::Vector3d>& verts, const std::vector<int>& tris) {
  overts = verts;
  otris = tris;
  object.BeginModel();
  AddToModel(object, verts, tris);
  object.EndModel();
}

void CollisionSystemPQP::GetCollisions(std::vector<unsigned int>& objectVertexToFace, std::vector<unsigned int>& edgeToEdge, std::vector<double>& edgeU, std::vector<Eigen::Vector3d>& moveEdge) {
  PQP_REAL translation[3];
  translation[0] = 0;
  translation[1] = 0;
  translation[2] = 0;
  PQP_REAL rotation[3][3];
  rotation[0][0] = 1;
  rotation[1][0] = 0;
  rotation[2][0] = 0;
  rotation[0][1] = 0;
  rotation[1][1] = 1;
  rotation[2][1] = 0;
  rotation[0][2] = 0;
  rotation[1][2] = 0;
  rotation[2][2] = 1;

  PQP_CollideResult cres;
  PQP_Collide(&cres, rotation, translation, &(ground), rotation, translation, &(object));
    int eToE = 0;
    int theirV = 0;
    int coPlane = 0;
  for (int j = 0; j < cres.NumPairs(); j++) {
    //determine possible vertex to face collision
    Eigen::Vector3d v0, v1, v2, u0, u1,u2;
    v0 = gverts[gtris[cres.pairs[j].id1 * 3]];
    v1 = gverts[gtris[cres.pairs[j].id1 * 3 + 1]];
    v2 = gverts[gtris[cres.pairs[j].id1 * 3 + 2]];
    u0 = overts[otris[cres.pairs[j].id2 * 3]];
    u1 = overts[otris[cres.pairs[j].id2 * 3 + 1]];
    u2 = overts[otris[cres.pairs[j].id2 * 3 + 2]];
    //First find line of plane intersection
    Eigen::Vector3d vn, un;
    vn = ((v1 - v0).cross(v2 - v1));
    vn.normalize();
    un = ((u1 - u0).cross(u2 - u1));
    un.normalize();
    double vp, up;
    vp = v0.dot(-1 * vn);
    up = u0.dot(-1 * un);
    Eigen::Vector3d dir;
    dir = vn.cross(un);
    if (dir.norm() < .0000001) {
      fprintf(stderr, "co planar brah\n");
      continue;
    }
    dir.normalize();
    Eigen::Matrix<double, 2, 3> m;
    m << vn[0], vn[1], vn[2],
        un[0], un[1], un[2];
    Eigen::Vector2d b;
    b << -1 * vp, -1* up;
    Eigen::Vector3d point = m.colPivHouseholderQr().solve(b);

    // get biggest component of line
    int bigComp = 0;
    if (fabs(dir[1]) > fabs(dir[0])) {
      bigComp = 1;
      if (fabs(dir[2]) > fabs(dir[1])) {
        bigComp = 2;
      }
    } else if (fabs(dir[2]) > fabs(dir[0])) {
      bigComp = 2;
    }
    // now get the edge collisions with line
    bool v0above = (v0 - u0).dot(un) > 0;
    bool v1above = (v1 - u0).dot(un) > 0;
    bool v2above = (v2 - u0).dot(un) > 0;
    Eigen::Vector3d ev1, ev2, ov;
    if (v0above == v1above && v2above != v0above) {
      ev1 = v0;
      ev2 = v1;
      ov = v2;
    } else if (v0above == v2above && v1above != v0above) {
      ev1 = v0;
      ev2 = v2;
      ov = v1;
    } else if (v2above == v1above && v2above != v0above) {
      ev1 = v2;
      ev2 = v1;
      ov = v0;
    } else {
      //fprintf(stderr, "just touching?\n");
      continue;
    }
    double evs1, evs2;
    evs1 = ((-1 * un).dot(ov - u0))/(un.dot(ev1 - ov));
    evs2 = ((-1 * un).dot(ov - u0))/(un.dot(ev2 - ov));


    double vi1 = ((evs1 * (ev1[bigComp] - ov[bigComp]) + ov[bigComp]) - point[bigComp])/ dir[bigComp];
    double vi2 = ((evs2 * (ev2[bigComp] - ov[bigComp]) + ov[bigComp]) - point[bigComp])/ dir[bigComp];
    
    bool u0above = (u0 - v0).dot(vn) > 0;
    bool u1above = (u1 - v0).dot(vn) > 0;
    bool u2above = (u2 - v0).dot(vn) > 0;
    Eigen::Vector3d eu1, eu2, ou;
    int utype;
    int utype2;
    if (u0above == u1above && u2above != u0above) {
      eu1 = u0;
      eu2 = u1;
      ou = u2;
      utype = 2;
      utype2 = 0;
      if (u0above == false) {
        //fprintf(stderr, "both under\n");
        utype += 3;
      }
    } else if (u0above == u2above && u1above != u0above) {
      eu1 = u0;
      eu2 = u2;
      ou = u1;
      utype = 1;
      utype2 = 0;
      if (u0above == false) {
        utype += 3;
        //fprintf(stderr, "both under\n");
      }
    } else if (u2above == u1above && u2above != u0above) {
      eu1 = u2;
      eu2 = u1;
      ou = u0;
      utype = 0;
      utype2 = 2;
      if (u2above == false) {
        utype += 3;
        //fprintf(stderr, "both under\n");
      }
    } else {
      //fprintf(stderr, "just touching?\n");
      continue;
    }
    double eus1, eus2;
    eus1 = ((-1 * vn).dot(ou - v0))/(vn.dot(eu1 - ou));
    eus2 = ((-1 * vn).dot(ou - v0))/(vn.dot(eu2 - ou));


    double ui1 = ((eus1 * (eu1[bigComp] - ou[bigComp]) + ou[bigComp]) - point[bigComp])/ dir[bigComp];
    double ui2 = ((eus2 * (eu2[bigComp] - ou[bigComp]) + ou[bigComp]) - point[bigComp])/ dir[bigComp];

    bool vflip = vi1 > vi2;
    bool uflip = ui1 > ui2;
    double temp;
    if (vflip) {
      temp = vi1;
      vi1 = vi2;
      vi2 = temp;
    }
    if (uflip) {
      temp = ui1;
      ui1 = ui2;
      ui2 = temp;
    }
    if (vi1 < ui1 && vi2 > ui2) {
      //fprintf(stderr, "Got one!\n");
      // correct vertexToFace!
      if (utype >= 3) {
        // the two other vertices are intersecting the face
        objectVertexToFace.push_back(cres.pairs[j].id2 * 3 + (utype +1)%3); // vertex
        objectVertexToFace.push_back(cres.pairs[j].id1 * 3); // first v of tri face
        objectVertexToFace.push_back(cres.pairs[j].id2 * 3 + (utype +2)%3); // vertex
        objectVertexToFace.push_back(cres.pairs[j].id1 * 3); // first v of tri face
      } else {
        objectVertexToFace.push_back(cres.pairs[j].id2 * 3 + utype); // vertex
        objectVertexToFace.push_back(cres.pairs[j].id1 * 3); // first v of tri face
      }

    } else {
      int oe = utype2;
      double esu;
      double movedir;
      if (vi1 > ui1) {
        movedir = vi1 - ui2;
      } else {
        movedir = vi2 - ui1;
      }
      if ((vi1 > ui1 && uflip) || (vi2 < ui2 && !uflip)) {
        oe = utype2;
        esu = eus1;
      } else if ((vi1 > ui1 && !uflip) || (vi2 < ui2 && uflip)) {
        oe = (utype2 + 1) %3;
        if (oe == utype%3) {
          oe = (oe + 1 )%3;
        }
        esu = eus2;
      } else {
        fprintf(stderr, "doesn't satisfy interval\n");
        continue;
      }
      edgeToEdge.push_back(cres.pairs[j].id2 * 3 + (utype %3));
      edgeToEdge.push_back(cres.pairs[j].id2 * 3 + oe);
      moveEdge.push_back(dir * movedir);
      edgeU.push_back(esu);
    }
      //fprintf(stderr, "not right vToF, vi1: %.3f vi2: %.3f ui1: %.3f ui2: %.3f\n", vi1, vi2, ui1, ui2);
      //fprintf(stderr, "dir[bigComp] = %.3f\n", dir[bigComp]);
    }

    //double v0[3], v1[3], v2[3];
    //double u0[3], u1[3], u2[3];
    //Eigen::Vector3d*p;
    // p = &(gverts[gtris[cres.pairs[j].id1 * 3]]);
    // v0[0] = (*p)[0];
    // v0[1] = (*p)[1];
    // v0[2] = (*p)[2];
    // p = &(gverts[gtris[cres.pairs[j].id1 * 3 + 1]]);
    // v1[0] = (*p)[0];
    // v1[1] = (*p)[1];
    // v1[2] = (*p)[2];
    // p = &(gverts[gtris[cres.pairs[j].id1 * 3 + 2]]);
    // v2[0] = (*p)[0];
    // v2[1] = (*p)[1];
    // v2[2] = (*p)[2];
    //p = &(overts[otris[cres.pairs[j].id2 * 3]]);
    // u0[0] = (*p)[0];
    // u0[1] = (*p)[1];
    // u0[2] = (*p)[2];
    //p = &(overts[otris[cres.pairs[j].id2 * 3 + 1]]);
    // u1[0] = (*p)[0];
    // u1[1] = (*p)[1];
    // u1[2] = (*p)[2];
    //p = &(overts[otris[cres.pairs[j].id2 * 3 + 2]]);
    // u2[0] = (*p)[0];
    // u2[1] = (*p)[1];
    // u2[2] = (*p)[2];
    //// ignoring coplanar triangles sorry
    ////

    //int ret;
    //if (ret = NoDivTriTriIsect(v0, v1, v2, u0, u1, u2)) {
    //  if (ret != 1) {
    //    fprintf(stderr, "got %i from ret, throwing out\n", ret);
    //    coPlane += 1;
    //  }
    //  // Figure out if it's edge to edge or vertex to edge

    //  if (vinter[0] > uinter[0] && vinter[1] < uinter[1]) {
    //    // vertex to face with u being the face
    //    // u is object so we don't handle it
    //    //fprintf(stderr, "Wrong vert to face\n");
    //    theirV += 1;
    //    continue;
    //  }
    //  if (uinter[0] > vinter[0] && uinter[1] < vinter[1]) {
    //    // vertex to face with v being the face
    //    // what is the vertex on u?
    //    //fprintf(stderr, "Correct vertex to face\n");
    //    objectVertexToFace.push_back(cres.pairs[j].id2 * 3 + tri_ue1[0]); // vertex
    //    objectVertexToFace.push_back(cres.pairs[j].id1 * 3); // first v of tri face
    //    continue;
    //  }
    //  eToE += 1;
    //  // prolly edge to edge
    //} else {
    //  fprintf(stderr, "didn't get collision!\n");
    //  if (TriContact(v0, v1, v2, u0, u1, u2)) {
    //    fprintf(stderr, "Their tri contact works :<\n");
    //  }
    //}
 // }
  if (cres.NumPairs() > 0) {
    //fprintf(stderr, "pairs: %i, eToE: %i, theirV: %i, coPlane: %i\n", cres.NumPairs(), eToE, theirV, coPlane);
  }
}

//http://fileadmin.cs.lth.se/cs/Personal/Tomas_Akenine-Moller/code/opttritri.txt
/* Triangle/triangle intersection test routine,
 * by Tomas Moller, 1997.
 * See article "A Fast Triangle-Triangle Intersection Test",
 * Journal of Graphics Tools, 2(2), 1997
 *
 * Updated June 1999: removed the divisions -- a little faster now!
 * Updated October 1999: added {} to CROSS and SUB macros 
 *
 * int NoDivTriTriIsect(float V0[3],float V1[3],float V2[3],
 *                      float U0[3],float U1[3],float U2[3])
 *
 * parameters: vertices of triangle 1: V0,V1,V2
 *             vertices of triangle 2: U0,U1,U2
 * result    : returns 1 if the triangles intersect, otherwise 0
 *
 */


/* if USE_EPSILON_TEST is true then we do a check:
         if |dv|<EPSILON then dv=0.0;
   else no check is done (which is less robust)
*/
//#define USE_EPSILON_TEST TRUE
#define EPSILON 0.000001


/* some macros */
#define CROSS(dest,v1,v2){                     \
              dest[0]=v1[1]*v2[2]-v1[2]*v2[1]; \
              dest[1]=v1[2]*v2[0]-v1[0]*v2[2]; \
              dest[2]=v1[0]*v2[1]-v1[1]*v2[0];}

#define DOT(v1,v2) (v1[0]*v2[0]+v1[1]*v2[1]+v1[2]*v2[2])

#define SUB(dest,v1,v2){         \
            dest[0]=v1[0]-v2[0]; \
            dest[1]=v1[1]-v2[1]; \
            dest[2]=v1[2]-v2[2];}

/* sort so that a<=b */
#define SORT(a,b,f)       \
             f = false; \
             if(a>b)    \
             {          \
               f = true; \
               double c; \
               c=a;     \
               a=b;     \
               b=c;     \
             }


/* this edge to edge test is based on Franlin Antonio's gem:
   "Faster Line Segment Intersection", in Graphics Gems III,
   pp. 199-202 */
#define EDGE_EDGE_TEST(V0,U0,U1)                      \
  Bx=U0[i0]-U1[i0];                                   \
  By=U0[i1]-U1[i1];                                   \
  Cx=V0[i0]-U0[i0];                                   \
  Cy=V0[i1]-U0[i1];                                   \
  f=Ay*Bx-Ax*By;                                      \
  d=By*Cx-Bx*Cy;                                      \
  if((f>0 && d>=0 && d<=f) || (f<0 && d<=0 && d>=f))  \
  {                                                   \
    e=Ax*Cy-Ay*Cx;                                    \
    if(f>0)                                           \
    {                                                 \
      if(e>=0 && e<=f) return 2;                      \
    }                                                 \
    else                                              \
    {                                                 \
      if(e<=0 && e>=f) return 2;                      \
    }                                                 \
  }

#define EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2) \
{                                              \
  double Ax,Ay,Bx,By,Cx,Cy,e,d,f;               \
  Ax=V1[i0]-V0[i0];                            \
  Ay=V1[i1]-V0[i1];                            \
  /* test edge U0,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U0,U1);                    \
  /* test edge U1,U2 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U1,U2);                    \
  /* test edge U2,U1 against V0,V1 */          \
  EDGE_EDGE_TEST(V0,U2,U0);                    \
}

#define POINT_IN_TRI(V0,U0,U1,U2)           \
{                                           \
  double a,b,c,d0,d1,d2;                     \
  /* is T1 completly inside T2? */          \
  /* check if V0 is inside tri(U0,U1,U2) */ \
  a=U1[i1]-U0[i1];                          \
  b=-(U1[i0]-U0[i0]);                       \
  c=-a*U0[i0]-b*U0[i1];                     \
  d0=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U2[i1]-U1[i1];                          \
  b=-(U2[i0]-U1[i0]);                       \
  c=-a*U1[i0]-b*U1[i1];                     \
  d1=a*V0[i0]+b*V0[i1]+c;                   \
                                            \
  a=U0[i1]-U2[i1];                          \
  b=-(U0[i0]-U2[i0]);                       \
  c=-a*U2[i0]-b*U2[i1];                     \
  d2=a*V0[i0]+b*V0[i1]+c;                   \
  if(d0*d1>0.0)                             \
  {                                         \
    if(d0*d2>0.0) return 2;                 \
  }                                         \
}

static int coplanar_tri_tri(double N[3],double V0[3],double V1[3],double V2[3],
                     double U0[3],double U1[3],double U2[3])
{
   double A[3];
   short i0,i1;
   /* first project onto an axis-aligned plane, that maximizes the area */
   /* of the triangles, compute indices: i0,i1. */
   A[0]=FABS(N[0]);
   A[1]=FABS(N[1]);
   A[2]=FABS(N[2]);
   if(A[0]>A[1])
   {
      if(A[0]>A[2])
      {
          i0=1;      /* A[0] is greatest */
          i1=2;
      }
      else
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
   }
   else   /* A[0]<=A[1] */
   {
      if(A[2]>A[1])
      {
          i0=0;      /* A[2] is greatest */
          i1=1;
      }
      else
      {
          i0=0;      /* A[1] is greatest */
          i1=2;
      }
    }

    /* test all edges of triangle 1 against the edges of triangle 2 */
    EDGE_AGAINST_TRI_EDGES(V0,V1,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V1,V2,U0,U1,U2);
    EDGE_AGAINST_TRI_EDGES(V2,V0,U0,U1,U2);

    /* finally, test if tri1 is totally contained in tri2 or vice versa */
    POINT_IN_TRI(V0,U0,U1,U2);
    POINT_IN_TRI(U0,V0,V1,V2);

    return 0;
}



#define NEWCOMPUTE_INTERVALS(VV0,VV1,VV2,D0,D1,D2,D0D1,D0D2,A,B,C,X0,X1, E1, E2) \
{ \
        if(D0D1>0.0f) \
        { \
                /* here we know that D0D2<=0.0 */ \
            /* that is D0, D1 are on the same side, D2 on the other or on the plane */ \
                A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
                E1[0] = 2; E1[1] = 0; E2[0] = 2; E2[1] = 1;\
        } \
        else if(D0D2>0.0f)\
        { \
                /* here we know that d0d1<=0.0 */ \
            A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
                E1[0] = 1; E1[1] = 0; E2[0] = 1; E2[1] = 2;\
        } \
        else if(D1*D2>0.0f || D0!=0.0f) \
        { \
                /* here we know that d0d1<=0.0 or that D0!=0.0 */ \
                A=VV0; B=(VV1-VV0)*D0; C=(VV2-VV0)*D0; X0=D0-D1; X1=D0-D2; \
                E1[0] = 0; E1[1] = 1; E2[0] = 0; E2[1] = 2;\
        } \
        else if(D1!=0.0f) \
        { \
                A=VV1; B=(VV0-VV1)*D1; C=(VV2-VV1)*D1; X0=D1-D0; X1=D1-D2; \
                E1[0] = 1; E1[1] = 0; E2[0] = 1; E2[1] = 2;\
        } \
        else if(D2!=0.0f) \
        { \
                A=VV2; B=(VV0-VV2)*D2; C=(VV1-VV2)*D2; X0=D2-D0; X1=D2-D1; \
                E1[0] = 2; E1[1] = 0; E2[0] = 2; E2[1] = 1;\
        } \
        else \
        { \
                /* triangles are coplanar */ \
                return coplanar_tri_tri(N1,V0,V1,V2,U0,U1,U2); \
        } \
}



static int NoDivTriTriIsect(double V0[3],double V1[3],double V2[3],
                     double U0[3],double U1[3],double U2[3])
{
  double E1[3],E2[3];
  double N1[3],N2[3],d1,d2;
  double du0,du1,du2,dv0,dv1,dv2;
  double D[3];
  double isect1[2], isect2[2];
  double du0du1,du0du2,dv0dv1,dv0dv2;
  short index;
  double vp0,vp1,vp2;
  double up0,up1,up2;
  double bb,cc,max;

  /* compute plane equation of triangle(V0,V1,V2) */
  SUB(E1,V1,V0);
  SUB(E2,V2,V0);
  CROSS(N1,E1,E2);
  d1=-DOT(N1,V0);
  /* plane equation 1: N1.X+d1=0 */

  /* put U0,U1,U2 into plane equation 1 to compute signed distances to the plane*/
  du0=DOT(N1,U0)+d1;
  du1=DOT(N1,U1)+d1;
  du2=DOT(N1,U2)+d1;

  /* coplanarity robustness check */
#if USE_EPSILON_TEST==TRUE
  if(FABS(du0)<EPSILON) du0=0.0;
  if(FABS(du1)<EPSILON) du1=0.0;
  if(FABS(du2)<EPSILON) du2=0.0;
#endif
  du0du1=du0*du1;
  du0du2=du0*du2;

  if(du0du1>0.0f && du0du2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute plane of triangle (U0,U1,U2) */
  SUB(E1,U1,U0);
  SUB(E2,U2,U0);
  CROSS(N2,E1,E2);
  d2=-DOT(N2,U0);
  /* plane equation 2: N2.X+d2=0 */

  /* put V0,V1,V2 into plane equation 2 */
  dv0=DOT(N2,V0)+d2;
  dv1=DOT(N2,V1)+d2;
  dv2=DOT(N2,V2)+d2;

#if USE_EPSILON_TEST==TRUE
  if(FABS(dv0)<EPSILON) dv0=0.0;
  if(FABS(dv1)<EPSILON) dv1=0.0;
  if(FABS(dv2)<EPSILON) dv2=0.0;
#endif

  dv0dv1=dv0*dv1;
  dv0dv2=dv0*dv2;

  if(dv0dv1>0.0f && dv0dv2>0.0f) /* same sign on all of them + not equal 0 ? */
    return 0;                    /* no intersection occurs */

  /* compute direction of intersection line */
  CROSS(D,N1,N2);

  /* compute and index to the largest component of D */
  max=(double)FABS(D[0]);
  index=0;
  bb=(double)FABS(D[1]);
  cc=(double)FABS(D[2]);
  if(bb>max) max=bb,index=1;
  if(cc>max) max=cc,index=2;

  /* this is the simplified projection onto L*/
  vp0=V0[index];
  vp1=V1[index];
  vp2=V2[index];

  up0=U0[index];
  up1=U1[index];
  up2=U2[index];

  /* compute interval for triangle 1 */
  double a,b,c,x0,x1;
  NEWCOMPUTE_INTERVALS(vp0,vp1,vp2,dv0,dv1,dv2,dv0dv1,dv0dv2,a,b,c,x0,x1, tri_ve1, tri_ve2);

  /* compute interval for triangle 2 */
  double d,e,f,y0,y1;
  NEWCOMPUTE_INTERVALS(up0,up1,up2,du0,du1,du2,du0du1,du0du2,d,e,f,y0,y1, tri_ue1, tri_ue2);

  double xx,yy,xxyy,tmp;
  xx=x0*x1;
  yy=y0*y1;
  xxyy=xx*yy;

  tmp=a*xxyy;
  isect1[0]=tmp+b*x1*yy;
  isect1[1]=tmp+c*x0*yy;

  tmp=d*xxyy;
  isect2[0]=tmp+e*xx*y1;
  isect2[1]=tmp+f*xx*y0;

  SORT(isect1[0],isect1[1], vflip);
  SORT(isect2[0],isect2[1], uflip);

  vinter[0] = isect1[0];
  vinter[1] = isect1[1];
  uinter[0] = isect2[0];
  uinter[1] = isect2[1];

  if(isect1[1]<isect2[0] || isect2[1]<isect1[0]) return 0;
  return 1;
}

static void VcrossV(PQP_REAL Vr[3], const PQP_REAL V1[3], const PQP_REAL V2[3])
{
  Vr[0] = V1[1]*V2[2] - V1[2]*V2[1];
  Vr[1] = V1[2]*V2[0] - V1[0]*V2[2];
  Vr[2] = V1[0]*V2[1] - V1[1]*V2[0];
}

static inline
PQP_REAL
Vlength(PQP_REAL V[3])
{
  return sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
}

static
inline
void
Vnormalize(PQP_REAL V[3])
{
  PQP_REAL d = (PQP_REAL)1.0 / sqrt(V[0]*V[0] + V[1]*V[1] + V[2]*V[2]);
  V[0] *= d;
  V[1] *= d;
  V[2] *= d;
}

static
inline
PQP_REAL
VdotV(const PQP_REAL V1[3], const PQP_REAL V2[3])
{
  return (V1[0]*V2[0] + V1[1]*V2[1] + V1[2]*V2[2]);
}

// TRIANGLE OVERLAP TEST
       
static
inline
PQP_REAL
max(PQP_REAL a, PQP_REAL b, PQP_REAL c)
{
  PQP_REAL t = a;
  if (b > t) t = b;
  if (c > t) t = c;
  return t;
}

static
inline
PQP_REAL
min(PQP_REAL a, PQP_REAL b, PQP_REAL c)
{
  PQP_REAL t = a;
  if (b < t) t = b;
  if (c < t) t = c;
  return t;
}

static
int
project6(PQP_REAL *ax, 
         PQP_REAL *p1, PQP_REAL *p2, PQP_REAL *p3, 
         PQP_REAL *q1, PQP_REAL *q2, PQP_REAL *q3)
{
  PQP_REAL P1 = VdotV(ax, p1);
  PQP_REAL P2 = VdotV(ax, p2);
  PQP_REAL P3 = VdotV(ax, p3);
  PQP_REAL Q1 = VdotV(ax, q1);
  PQP_REAL Q2 = VdotV(ax, q2);
  PQP_REAL Q3 = VdotV(ax, q3);
  
  PQP_REAL mx1 = max(P1, P2, P3);
  PQP_REAL mn1 = min(P1, P2, P3);
  PQP_REAL mx2 = max(Q1, Q2, Q3);
  PQP_REAL mn2 = min(Q1, Q2, Q3);

  if (mn1 > mx2) return 0;
  if (mn2 > mx1) return 0;
  return 1;
}
static int 
TriContact(PQP_REAL *P1, PQP_REAL *P2, PQP_REAL *P3,
           PQP_REAL *Q1, PQP_REAL *Q2, PQP_REAL *Q3) 
{

  // One triangle is (p1,p2,p3).  Other is (q1,q2,q3).
  // Edges are (e1,e2,e3) and (f1,f2,f3).
  // Normals are n1 and m1
  // Outwards are (g1,g2,g3) and (h1,h2,h3).
  //  
  // We assume that the triangle vertices are in the same coordinate system.
  //
  // First thing we do is establish a new c.s. so that p1 is at (0,0,0).

  PQP_REAL p1[3], p2[3], p3[3];
  PQP_REAL q1[3], q2[3], q3[3];
  PQP_REAL e1[3], e2[3], e3[3];
  PQP_REAL f1[3], f2[3], f3[3];
  PQP_REAL g1[3], g2[3], g3[3];
  PQP_REAL h1[3], h2[3], h3[3];
  PQP_REAL n1[3], m1[3];

  PQP_REAL ef11[3], ef12[3], ef13[3];
  PQP_REAL ef21[3], ef22[3], ef23[3];
  PQP_REAL ef31[3], ef32[3], ef33[3];
  
  p1[0] = P1[0] - P1[0];  p1[1] = P1[1] - P1[1];  p1[2] = P1[2] - P1[2];
  p2[0] = P2[0] - P1[0];  p2[1] = P2[1] - P1[1];  p2[2] = P2[2] - P1[2];
  p3[0] = P3[0] - P1[0];  p3[1] = P3[1] - P1[1];  p3[2] = P3[2] - P1[2];
  
  q1[0] = Q1[0] - P1[0];  q1[1] = Q1[1] - P1[1];  q1[2] = Q1[2] - P1[2];
  q2[0] = Q2[0] - P1[0];  q2[1] = Q2[1] - P1[1];  q2[2] = Q2[2] - P1[2];
  q3[0] = Q3[0] - P1[0];  q3[1] = Q3[1] - P1[1];  q3[2] = Q3[2] - P1[2];
  
  e1[0] = p2[0] - p1[0];  e1[1] = p2[1] - p1[1];  e1[2] = p2[2] - p1[2];
  e2[0] = p3[0] - p2[0];  e2[1] = p3[1] - p2[1];  e2[2] = p3[2] - p2[2];
  e3[0] = p1[0] - p3[0];  e3[1] = p1[1] - p3[1];  e3[2] = p1[2] - p3[2];

  f1[0] = q2[0] - q1[0];  f1[1] = q2[1] - q1[1];  f1[2] = q2[2] - q1[2];
  f2[0] = q3[0] - q2[0];  f2[1] = q3[1] - q2[1];  f2[2] = q3[2] - q2[2];
  f3[0] = q1[0] - q3[0];  f3[1] = q1[1] - q3[1];  f3[2] = q1[2] - q3[2];
  
  VcrossV(n1, e1, e2);
  VcrossV(m1, f1, f2);

  VcrossV(g1, e1, n1);
  VcrossV(g2, e2, n1);
  VcrossV(g3, e3, n1);
  VcrossV(h1, f1, m1);
  VcrossV(h2, f2, m1);
  VcrossV(h3, f3, m1);

  VcrossV(ef11, e1, f1);
  VcrossV(ef12, e1, f2);
  VcrossV(ef13, e1, f3);
  VcrossV(ef21, e2, f1);
  VcrossV(ef22, e2, f2);
  VcrossV(ef23, e2, f3);
  VcrossV(ef31, e3, f1);
  VcrossV(ef32, e3, f2);
  VcrossV(ef33, e3, f3);
  
  // now begin the series of tests

  if (!project6(n1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(m1, p1, p2, p3, q1, q2, q3)) return 0;
  
  if (!project6(ef11, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef12, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef13, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef21, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef22, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef23, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef31, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef32, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(ef33, p1, p2, p3, q1, q2, q3)) return 0;

  if (!project6(g1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(g2, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(g3, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(h1, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(h2, p1, p2, p3, q1, q2, q3)) return 0;
  if (!project6(h3, p1, p2, p3, q1, q2, q3)) return 0;

  return 1;
}
