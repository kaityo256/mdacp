#include <iostream>
#include "mpistream.h"
#include "simulationinfo.h"
//----------------------------------------------------------------------
SimulationInfo::SimulationInfo(Parameter &param, int *grid_size) {
  BufferLength = 0.3;
  SearchLength = CUTOFF_LENGTH + BufferLength;
  BaseDir = ".";
  TimeStep = param.GetDoubleDef("TimeStep", 0.001);
  IsPeriodic = param.GetBooleanDef("IsPeriodic", true);
  AimedTemperature = param.GetDoubleDef("AimedTemperature", 1.0);
  ControlTemperature = param.GetBooleanDef("ControlTemperature", false);
  HeatbathTau = param.GetDoubleDef("HeatbathTau", 0.1);
  HeatbathGamma = param.GetDoubleDef("HeatbathGamma", 0.1);
  BaseDir = param.GetStringDef("BaseDir", ".");
  SortParticle = param.GetBooleanDef("SortParticle", false);

  std::string hbtype = param.GetStringDef("HeatbathType", "NoseHoover");

  if (hbtype == "NoseHoover") {
    HeatbathType = HT_NOSEHOOVER;
  } else if (hbtype == "Langevin") {
    HeatbathType = HT_LANGEVIN;
  } else {
    mout << "Error: Unknown Heatbath type " << hbtype << std::endl;
  }

  //For PairList
  CheckListLength = param.GetIntegerDef("CheckListLength", 500);

  if (param.Contains("SystemSize")) {
    double ss = param.GetDouble("SystemSize");
    L[X] = ss;
    L[Y] = ss;
    L[Z] = ss;
  } else if (param.Contains("UnitLength")) {
    const double ul = param.GetDouble("UnitLength");
    L[X] = ul * grid_size[X];
    L[Y] = ul * grid_size[Y];
    L[Z] = ul * grid_size[Z];
  }
  if (param.Contains("SystemSizeX")) {
    L[X] = param.GetDouble("SystemSizeX");
  }
  if (param.Contains("SystemSizeY")) {
    L[Y] = param.GetDouble("SystemSizeY");
  }
  if (param.Contains("SystemSizeZ")) {
    L[Z] = param.GetDouble("SystemSizeZ");
  }
}
//----------------------------------------------------------------------
void
SimulationInfo::ShowAll(unsigned long int pn) {
  double density = static_cast<double>(pn) / (L[X] * L[Y] * L[Z]);
  mout << "# Number of Particles = " << pn << std::endl;
  mout << "# System Size = (" << L[X] << ",";
  mout << L[Y] << "," << L[Z] << ")" << std::endl;
  mout << "# Density = " << density << std::endl;
  mout << "# Cutoff Length  = " << CUTOFF_LENGTH << std::endl;
  mout << "# TimeStep = " << TimeStep << std::endl;
  mout << "# ControlTemperature = " << (ControlTemperature ? "yes" : "no") << std::endl;
  if (ControlTemperature) {
    mout << "# AimedTemperature =" << AimedTemperature << std::endl;
    if (HT_NOSEHOOVER == HeatbathType) {
      mout << "# HeatbathType = NoseHoover" << std::endl;
      mout << "# HeatbathTau = " << HeatbathTau << std::endl;
    } else if (HT_LANGEVIN == HeatbathType) {
      mout << "# HeatbathType = Langevin" << std::endl;
      mout << "# HeatbathGamma = " << HeatbathGamma << std::endl;
    }
  }
  mout << "# IsPeriodic = " << (IsPeriodic ? "yes" : "no") << std::endl;
}
//----------------------------------------------------------------------
void
SimulationInfo::AdjustPeriodicBoundary(double &dx, double &dy, double &dz) {
  if (dx < -L[X] * 0.5) {
    dx += L[X];
  } else if (dx > L[X] * 0.5) {
    dx -= L[X];
  }
  if (dy < -L[Y] * 0.5) {
    dy += L[Y];
  } else if (dy > L[Y] * 0.5) {
    dy -= L[Y];
  }
  if (dz < -L[Z] * 0.5) {
    dz += L[Z];
  } else if (dz > L[Z] * 0.5) {
    dz -= L[Z];
  }
}
//---------------------------------------------------------------------
