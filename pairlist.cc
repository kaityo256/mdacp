//----------------------------------------------------------------------
// Check Validity of Pair-list
//----------------------------------------------------------------------
#include <iostream>
#include <vector>
#include "mpistream.h"
#include "communicator.h"
#include "pairlist.h"
//----------------------------------------------------------------------
void
PairList::Init(Variables *vars, SimulationInfo *sinfo) {
  initialized = true;
  is_fresh = true;
  buffer_length = sinfo->BufferLength;
  lifetime = 0;
  double (*q)[D] = vars->q;
  const int pn = vars->GetTotalParticleNumber();
  for (int i = 0; i < pn; i++) {
    qb_old[i][X] = q[i][X];
    qb_old[i][Y] = q[i][Y];
    qb_old[i][Z] = q[i][Z];
  }
}
//----------------------------------------------------------------------
bool
PairList::IsPairListExpired(Variables *vars, SimulationInfo *sinfo) {
  if (!initialized) {
    return true;
  }
  is_fresh = false;
  double (*p)[D] = vars->p;

  lifetime++;
  if (buffer_length > 0.0) {
    UpdatePairListValidity(vars->GetParticleNumber(), p, sinfo);
  }

  //return (buffer_length < 0.0);

  bool expired = false;
  if (buffer_length < 0.0) {
    expired = CheckByDisplacement(vars, sinfo);
  }
  return expired;
}
//----------------------------------------------------------------------
void
PairList::UpdatePairListValidity(const int pn, double p[N][D], SimulationInfo *sinfo) {
  double max_velocity = p[0][X] * p[0][X] + p[0][Y] * p[0][Y] + p[0][Z] * p[0][Z];
  int index = 0;
  for (int i = 1; i < pn; i++) {
    double v = p[i][X] * p[i][X] + p[i][Y] * p[i][Y] + p[i][Z] * p[i][Z];
    if (v > max_velocity) {
      max_velocity = v;
      index = i;
    }
  }
  max_velocity = sqrt(max_velocity);
  if (max_velocity * 2.0 * sinfo->TimeStep > sinfo->BufferLength) {
    show_error("Too fast particles exists. Try smaller time step");
    printf("%d :max_velocity = %f\n", index, max_velocity);
    exit(1);
  }
  buffer_length = buffer_length - max_velocity * 2.0 * sinfo->TimeStep;
}
//----------------------------------------------------------------------
bool
PairList::CheckByDisplacement(Variables *vars, SimulationInfo *sinfo) {
  const int pn = vars->GetTotalParticleNumber();
  double (*q)[D] = vars->q;

  double dr2_max = 0.0;
  for (int i = 0; i < pn; i++) {
    double dx = qb_old[i][X] - q[i][X];
    double dy = qb_old[i][Y] - q[i][Y];
    double dz = qb_old[i][Z] - q[i][Z];
    double dr2 = dx * dx + dy * dy + dz * dz;
    disp[i] = dr2;
    if (dr2_max < dr2) {
      dr2_max = dr2;
    }
  }


  const int Nl = (pn > sinfo->CheckListLength) ? sinfo->CheckListLength : pn;
  std::vector<int> v;

  const double c_distance = sinfo->BufferLength - sqrt(dr2_max);
  const double c_d2 = c_distance * c_distance;
  for (int i = 0; i < pn; i++) {
    if (disp[i] > c_d2) {
      v.push_back(i);
      if (static_cast<int>(v.size()) >= Nl) {
        return true;
      }
    }
  }

  const double SL2 = sinfo->SearchLength * sinfo->SearchLength;
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;

  const int l = v.size();

  for (int i = 0; i < l - 1; i++) {
    for (int j = i + 1; j < l; j++) {
      int i1 = v[i];
      int i2 = v[j];
      double dx = qb_old[i1][X] - qb_old[i2][X];
      double dy = qb_old[i1][Y] - qb_old[i2][Y];
      double dz = qb_old[i1][Z] - qb_old[i2][Z];
      double dr2 = dx * dx + dy * dy + dz * dz;
      if (dr2 < SL2) continue;
      dx = q[i1][X] - q[i2][X];
      dy = q[i1][Y] - q[i2][Y];
      dz = q[i1][Z] - q[i2][Z];
      dr2 = dx * dx + dy * dy + dz * dz;
      if (dr2 < CL2) return true;
    }
  }
  return false;
}
//----------------------------------------------------------------------
