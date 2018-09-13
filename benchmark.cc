//----------------------------------------------------------------------
#include "benchmark.h"
#include "mpistream.h"
#include "communicator.h"
#include "confmaker.h"
#ifdef FX10
#include <fj_tool/fipp.h>
#endif
//----------------------------------------------------------------------
Benchmark bench;
//----------------------------------------------------------------------
void
Benchmark::Run(MDManager *mdm) {
  Parameter *param = mdm->GetParameter();
  //SimpleConfigurationMaker c(param);
  ExtractConfigurationMaker c(param);
  mdm->ExecuteAll(&c);
  const double v0 = param->GetDoubleDef("InitialVelocity", 1.0);
  mdm->SetInitialVelocity(v0);
  const int T_LOOP = param->GetIntegerDef("ThermalizeLoop", 150);
  const int LOOP = param->GetIntegerDef("TotalLoop", 1000);
  const int OBSERVE_LOOP = param->GetIntegerDef("ObserveLoop", 100);
  mdm->ShowSystemInformation();
  mdm->MakePairList();
  for (int i = 0; i < T_LOOP; i++) {
    mdm->Calculate();
  }
  double start_time = Communicator::GetTime();
#ifdef FX10
  fipp_start();
#endif
  for (int i = 0; i < LOOP; i++) {
    mdm->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mout << mdm->GetSimulationTime();
      mout << " " << mdm->Temperature();
      mout << " " << mdm->Pressure();
      mout << " " << mdm->TotalEnergy();
      mout << " #observe" << std::endl;
    }
  }
#ifdef FX10
  fipp_stop();
#endif
  double sec = Communicator::GetTime() - start_time;
  const unsigned long int pn = mdm->GetTotalParticleNumber();
  double mups = static_cast<double>(LOOP);
  mups = mups * static_cast<double>(pn) / sec / 1.0e6;
  mout << "# N = " << pn << " ";
  mout << sec << " [SEC] ";
  mout << mups << " [MUPS]" << std::endl;
}
//----------------------------------------------------------------------
