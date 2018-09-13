//----------------------------------------------------------------------
#ifndef variables_h
#define variables_h
//----------------------------------------------------------------------
#include <iostream>
#include "mdconfig.h"
#include "simulationinfo.h"
//----------------------------------------------------------------------
class Variables {
private:
  int particle_number;
  int total_particle_number;
  double C0, C2;
public:
  Variables(void);
  int type[N];
  __attribute__((aligned(64))) double q[N][D];
  __attribute__((aligned(64))) double p[N][D];
  double Zeta;
  double SimulationTime;
  double GetC0(void) {return C0;}
  double GetC2(void) {return C2;}
  int GetParticleNumber(void) {return particle_number;}
  int GetTotalParticleNumber(void) {return total_particle_number;}

  void AddParticle(double x[D], double v[D], int t = 0);
  void SaveToStream(std::ostream &fs);
  void LoadFromStream(std::istream &fs);
  void SaveConfiguration(std::ostream &fs);
  void SetParticleNumber(int pn) {particle_number = pn;}
  void SetTotalParticleNumber(int pn) {total_particle_number = pn;}
  void AdjustPeriodicBoundary(SimulationInfo *sinfo);
  void SetInitialVelocity(double V0, const int id);

  void ChangeScale(double alpha);

  // For Debug
  void Shuffle(void);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
