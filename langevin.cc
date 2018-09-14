//----------------------------------------------------------------------
#include "langevin.h"
#include "mpistream.h"
#include "communicator.h"
#include "confmaker.h"
//----------------------------------------------------------------------
Langevin langevin;
//----------------------------------------------------------------------
void
Langevin::Run(MDManager *mdm) {
  Parameter *param = mdm->GetParameter();
  ConfigurationMaker c(param);
  mdm->ExecuteAll(&c);
  const double v0 = param->GetDoubleDef("InitialVelocity", 1.0);
  mdm->SetInitialVelocity(v0);
  const int LOOP = param->GetIntegerDef("TotalLoop", 1000);
  const int OBSERVE_LOOP = param->GetIntegerDef("ObserveLoop", 100);
  const int seed = param->GetIntegerDef("Seed",1);
  mdm->SetSeed(seed);
  mdm->SetControlTemperature(true);
  mdm->ShowSystemInformation();
  mdm->MakePairList();
  for (int i = 0; i < LOOP; i++) {
    mdm->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mout << mdm->GetSimulationTime();
      mout << " " << mdm->Temperature();
      mout << " " << mdm->Pressure();
      mout << " " << mdm->TotalEnergy() << std::endl;
    }
  }
}
//----------------------------------------------------------------------
