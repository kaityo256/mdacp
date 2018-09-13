//----------------------------------------------------------------------
#ifndef observe_h
#define observe_h
//----------------------------------------------------------------------
#include "variables.h"
#include "meshlist.h"
#include "simulationinfo.h"
class IntegerObserver {
public:
  virtual int Observe(Variables * vars, MeshList *mesh) = 0;
};
//----------------------------------------------------------------------
class DoubleObserver {
public:
  virtual double Observe(Variables * vars, MeshList *mesh) = 0;
};
//----------------------------------------------------------------------
class ParticleNumberObserver : public IntegerObserver {
public:
  int Observe(Variables *vars, MeshList *mesh = NULL) {return vars->GetParticleNumber();}
};
//----------------------------------------------------------------------
class KineticEnergyObserver : public DoubleObserver {
public:
  double Observe(Variables * vars, MeshList *mesh);
};
//----------------------------------------------------------------------
class PotentialEnergyObserver : public DoubleObserver {
private:
  SimulationInfo *sinfo;
public:
  PotentialEnergyObserver(SimulationInfo *sinfo_) {sinfo = sinfo_;}
  double Observe(Variables * vars, MeshList *mesh);
};
//----------------------------------------------------------------------
class VirialObserver : public DoubleObserver {
public:
  double Observe(Variables * vars, MeshList *mesh);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
