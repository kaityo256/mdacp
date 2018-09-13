//----------------------------------------------------------------------
// Check Validity of Pair-list
//----------------------------------------------------------------------
#ifndef pairlist_h
#define pairlist_h
#include "mdconfig.h"
#include "meshlist.h"
#include "variables.h"

class PairList {
private:
  int lifetime;
  bool initialized;
  bool is_fresh;
  double buffer_length;
  double qb_old[N][D];
  double disp[N];

  bool CheckByDisplacement(Variables *vars, SimulationInfo *sinfo);
  void UpdatePairListValidity(const int pn, double p[N][D], SimulationInfo *sinfo);
  static bool Compare(const int &i, const int &j);

public:

  PairList(void) {initialized = false;}
  void Init(Variables *vars, SimulationInfo *sinfo);
  bool IsFresh(void) {return is_fresh;}
  void SetFresh(bool b) {is_fresh = b;}
  bool IsPairListExpired(Variables *vars, SimulationInfo *sinfo);
  int GetLifeTime(void) {return lifetime;}
};
//----------------------------------------------------------------------
#endif
