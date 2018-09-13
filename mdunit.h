//----------------------------------------------------------------------
#ifndef mdunit_h
#define mdunit_h
#include <stdio.h>
#include <vector>
#include <random>
#include "mdconfig.h"
#include "parainfo.h"
#include "simulationinfo.h"
#include "mdrect.h"
#include "parameter.h"
#include "variables.h"
#include "meshlist.h"
#include "pairlist.h"
#include "observer.h"
#include "fcalculator.h"
//----------------------------------------------------------------------
class MDUnit;
//----------------------------------------------------------------------
class Executor {
public:
  virtual void Execute(MDUnit *) {};
};
//----------------------------------------------------------------------
class MDUnit {
private:
  Variables *vars;
  SimulationInfo *sinfo;
  ParaInfo *pinfo;
  MeshList *mesh;
  PairList *plist;
  const int id;
  std::vector<int> border_particles[MAX_DIR];
  MDRect myrect;
  std::mt19937 mt;
public:
  MDUnit (int id_, SimulationInfo *si, ParaInfo *pi);
  ~MDUnit(void);
  std::vector<ParticleInfo> send_buffer;
  int GetID(void) {return id;}
  MDRect * GetRect(void) {return &myrect;}
  Variables *GetVariables(void) {return vars;}
  void SaveConfiguration(void);
  void SaveAsCdview(std::ofstream &ofs);
  void AddParticle(double x[D], double v[D], int type = 0);
  void AddParticle(double x[D], int type = 0);
  double * GetSystemSize(void) {return sinfo->L;}
  double GetTimeStep(void) {return sinfo->TimeStep;}
  void SetInitialVelocity(double v0) {vars->SetInitialVelocity(v0, id);}
  int GetParticleNumber(void) {return vars->GetParticleNumber();}
  int GetTotalParticleNumber(void) {return vars->GetTotalParticleNumber();}
  void SetTotalParticleNumber(int n) {vars->SetTotalParticleNumber(n);}

  void CalculateForce(void) {ForceCalculator::CalculateForce(vars, mesh, sinfo);}
  void UpdatePositionHalf(void) {ForceCalculator::UpdatePositionHalf(vars, sinfo);}
  void HeatbathZeta(double t) {ForceCalculator::HeatbathZeta(vars, t, sinfo);}
  void HeatbathMomenta(void) {ForceCalculator::HeatbathMomenta(vars, sinfo);}
  void SetSeed(const int seed){mt.seed(seed + GetID());}
  void Langevin(void) {ForceCalculator::Langevin(vars, sinfo, mt);}

  void MakeBufferForSendingParticle(const int dir);
  void FindBorderParticles(const int dir);
  void MakeBufferForBorderParticles(const int dir);
  void ReceiveParticles(std::vector<ParticleInfo> &recv_buffer);
  void AdjustPeriodicBoundary(void) {vars->AdjustPeriodicBoundary(sinfo);}
  void ReceiveBorderParticles(std::vector<ParticleInfo> &recv_buffer);
  double ObserveDouble(DoubleObserver *obs) {return obs->Observe(vars, mesh);}
  int IntegerDouble(IntegerObserver *obs) {return obs->Observe(vars, mesh);}
  void Execute(Executor *ex) {ex->Execute(this);}
  void MakePairList(void);
  void ShowPairs(void) {mesh->ShowPairs();}
  bool IsPairListExpired(void) {return plist->IsPairListExpired(vars, sinfo);}
  //
  void ChangeScale(double alpha);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
