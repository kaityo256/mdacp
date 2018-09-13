//----------------------------------------------------------------------
#ifndef fcalculator_h
#define fcalculator_h
//----------------------------------------------------------------------
#include <random>
#include "mdconfig.h"
#include "variables.h"
#include "meshlist.h"
#include "simulationinfo.h"
//----------------------------------------------------------------------
namespace ForceCalculator {
void CalculateForceBruteforce(Variables *vars, SimulationInfo *sinfo);
void CalculateForceUnroll(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
void CalculateForceSorted(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
void CalculateForceNext(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
void CalculateForcePair(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
void CalculateForceReactless(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);

void CalculateForceReactlessSIMD(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
void CalculateForceReactlessSIMD_errsafe(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
#ifdef AVX2
void CalculateForceAVX2(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
#endif
void UpdatePositionHalf(Variables *vars, SimulationInfo *sinfo);
void CalculateForce(Variables *vars, MeshList *mesh, SimulationInfo *sinfo);
void HeatbathZeta(Variables *vars, double ct, SimulationInfo *sinfo);
void HeatbathMomenta(Variables *vars, SimulationInfo *sinfo);
void Langevin(Variables *vars, SimulationInfo *sinfo, std::mt19937 &mt);
};
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------
