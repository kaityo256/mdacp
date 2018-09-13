//---------------------------------------------------------------------
#ifndef simulationinfo_h
#define simulationinfo_h
//---------------------------------------------------------------------
#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "parameter.h"
#include "mdconfig.h"
//---------------------------------------------------------------------
class SimulationInfo {
private:

public:
  SimulationInfo(Parameter &param, int grid_size[D]);
  void ShowAll(unsigned long int pn);
  void AdjustPeriodicBoundary(double &x, double &y, double &z);

  //Public Parameters
  double L[D];
  double TimeStep;
  double HeatbathTau;
  double HeatbathGamma;
  double AimedTemperature;
  bool ControlTemperature;
  bool SortParticle;
  std::string BaseDir;
  double SearchLength;
  double BufferLength;
  bool IsPeriodic;
  int CheckListLength;
  int HeatbathType;

};
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------
