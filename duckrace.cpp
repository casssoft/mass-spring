#include "duckrace.h"
#include "particle_system.h"

#define P_SIZE m->particles.size()
namespace {

ParticleSystem* m;
void addleg(int from, double tx, double ty, double iMass) {
  m->particles.emplace_back();
  Particle*p = &(m->particles[m->particles.size()-1]);

  p->x << tx*2, ty*2, 0;
  p->iMass = iMass;
  p->v << 0.0, 0.0, 0.0;
}

void addspring(int from, int to, double L, double k) {
  m->springs.emplace_back();
  Spring*s = &(m->springs[m->springs.size()-1]);

  s->from = from;
  s->to = to;
  s->L = L*2;
  s->k = k*9000;
  s->c = (k*9000)/50;
}
};
void Duckrace::MakeBall(ParticleSystem* ps, double x1,double y1) {
  m = ps;
	addleg(-1,x1,y1,1);//0
	addleg(P_SIZE+7,x1,y1+10,1);//1

	addleg(P_SIZE-1,x1+7,y1+7,1);//2
	addleg(P_SIZE-1,x1+10,y1,1);//3
	addleg(P_SIZE-1,x1+7,y1-7,1);//4
	addleg(P_SIZE-1,x1,y1-10,1);//5
	addleg(P_SIZE-1,x1-7,y1-7,1);//6
	addleg(P_SIZE-1,x1-10,y1,1);//7
	addleg(P_SIZE-1,x1-7,y1+7,1);//8
	
	addspring(P_SIZE-9,P_SIZE-8,10,.008);
	addspring(P_SIZE-9,P_SIZE-7,10,.008);
	addspring(P_SIZE-9,P_SIZE-6,10,.008);
	addspring(P_SIZE-9,P_SIZE-5,10,.008);
	addspring(P_SIZE-9,P_SIZE-4,10,.008);	
	addspring(P_SIZE-9,P_SIZE-3,10,.008);
	addspring(P_SIZE-9,P_SIZE-2,10,.008);
	addspring(P_SIZE-9,P_SIZE-1,10,.008);
	
	
	addspring(P_SIZE-8,P_SIZE-7,7.6,.01);
	addspring(P_SIZE-7,P_SIZE-6,7.6,.01);
	addspring(P_SIZE-6,P_SIZE-5,7.6,.01);
	addspring(P_SIZE-5,P_SIZE-4,7.6,.01);
	addspring(P_SIZE-4,P_SIZE-3,7.6,.01);	
	addspring(P_SIZE-3,P_SIZE-2,7.6,.01);
	addspring(P_SIZE-2,P_SIZE-1,7.6,.01);
	addspring(P_SIZE-8,P_SIZE-1,7.6,.01);

	addspring(P_SIZE-8,P_SIZE-6,14,.008);
	addspring(P_SIZE-7,P_SIZE-5,14,.008);	
	addspring(P_SIZE-6,P_SIZE-4,14,.008);	
	addspring(P_SIZE-5,P_SIZE-3,14,.008);	
	addspring(P_SIZE-4,P_SIZE-2,14,.008);
	addspring(P_SIZE-3,P_SIZE-1,14,.008);
	addspring(P_SIZE-2,P_SIZE-8,14,.008);
	addspring(P_SIZE-1,P_SIZE-7,14,.008);

	addspring(P_SIZE-8,P_SIZE-4,20,.013);
	addspring(P_SIZE-7,P_SIZE-3,20,.013);
	addspring(P_SIZE-6,P_SIZE-2,20,.013);
	addspring(P_SIZE-5,P_SIZE-1,20,.013);
	
}
