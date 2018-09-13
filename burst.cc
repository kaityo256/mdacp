//----------------------------------------------------------------------
#include <fstream>
#include "burst.h"
#include "mpistream.h"
#include "communicator.h"
#include "confmaker.h"
//----------------------------------------------------------------------
Burst burst;
//----------------------------------------------------------------------
class DropletMaker : public Executor {
private:
  Parameter *param;
public:
  DropletMaker(Parameter *param_) {
    param = param_;
  };
  void AddSub(double x[D], MDUnit *mdu) {
    double *L = mdu->GetSystemSize();
    const double c = L[X] * 0.5;
    const double radius = param->GetDoubleDef("Radius", 10);;
    const double dx = x[X] - c;
    const double dy = x[Y] - c;
    const double dz = x[Z] - c;
    const double r2 = dx * dx + dy * dy + dz * dz;
    if (r2 < radius * radius) {
      mdu->AddParticle(x);
    }
  };
  void Execute(MDUnit *mdu) {
    double *L = mdu->GetSystemSize();
    const double density = param->GetDoubleDef("Density", 0.5);;
    const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
    const double hs = s * 0.5;
    int sx = static_cast<int>(L[X] / s);
    int sy = static_cast<int>(L[Y] / s);
    int sz = static_cast<int>(L[Z] / s);
    double x[D];
    const double e = 0.0000001;
    for (int iz = 0; iz < sz; iz++) {
      for (int iy = 0; iy < sy; iy++) {
        for (int ix = 0; ix < sx; ix++) {
          x[X] = static_cast<double>(ix) * s + e;
          x[Y] = static_cast<double>(iy) * s + e;
          x[Z] = static_cast<double>(iz) * s + e;
          AddSub(x, mdu);

          x[X] = static_cast<double>(ix) * s + e;
          x[Y] = static_cast<double>(iy) * s + hs + e;
          x[Z] = static_cast<double>(iz) * s + hs + e;
          AddSub(x, mdu);

          x[X] = static_cast<double>(ix) * s + hs + e;
          x[Y] = static_cast<double>(iy) * s + e;
          x[Z] = static_cast<double>(iz) * s + hs + e;
          AddSub(x, mdu);

          x[X] = static_cast<double>(ix) * s + hs + e;
          x[Y] = static_cast<double>(iy) * s + hs + e;
          x[Z] = static_cast<double>(iz) * s + e;
          AddSub(x, mdu);
        }
      }
    }
  };
};
//----------------------------------------------------------------------
class ShellForce : public Executor {
private:
  Parameter *param;
public:
  ShellForce(Parameter *param_) {
    param = param_;
  };
  void Execute(MDUnit *mdu) {
    double *L = mdu->GetSystemSize();
    const double c = L[X] * 0.5;
    const double radius = param->GetDoubleDef("Radius", 10) + 1.0;
    Variables *vars = mdu->GetVariables();
    const int pn = mdu->GetVariables()->GetParticleNumber();
    double (*q)[D] = mdu->GetVariables()->q;
    double (*p)[D] = mdu->GetVariables()->p;
    const double C2 = vars->GetC2();
    const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
    const double dt = 0.005;
    for (int i = 0; i < pn; i++) {
      const double dx = q[i][X] - c;
      const double dy = q[i][Y] - c;
      const double dz = q[i][Z] - c;
      const double r = radius - sqrt(dx * dx + dy * dy + dz * dz);
      const double r2 = r * r;
      if ( r2 > CL2) {
        continue;
      }
      const double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      p[i][X] += df * dx;
      p[i][Y] += df * dy;
      p[i][Z] += df * dz;
    }
  };
};
//----------------------------------------------------------------------
void
Burst::Run(MDManager *mdm) {
  Parameter *param = mdm->GetParameter();
  DropletMaker c(param);
  mdm->ExecuteAll(&c);
  const double v0 = param->GetDoubleDef("InitialVelocity", 1.0);
  mdm->SetInitialVelocity(v0);
  const int T_LOOP = param->GetIntegerDef("ThermalizeLoop", 150);
  const int LOOP = param->GetIntegerDef("TotalLoop", 1000);
  const int OBSERVE_LOOP = param->GetIntegerDef("ObserveLoop", 100);
  mdm->ShowSystemInformation();
  mdm->MakePairList();
  mdm->SetControlTemperature(true);
  for (int i = 0; i < T_LOOP; i++) {
    mdm->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mdm->SaveAsCdviewSequential();
    }
  }
  mdm->SetControlTemperature(false);
  double start_time = Communicator::GetTime();
  for (int i = 0; i < LOOP; i++) {
    mdm->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mout << mdm->GetSimulationTime();
      mout << " " << mdm->Temperature();
      mout << " " << mdm->Pressure();
      mout << " " << mdm->TotalEnergy();
      mout << "# observe" << std::endl;
      mdm->SaveAsCdviewSequential();
    }
  }
  double sec = Communicator::GetTime() - start_time;
  const unsigned long int pn = mdm->GetTotalParticleNumber();
  double mups = static_cast<double>(LOOP);
  mups = mups * static_cast<double>(pn) / sec / 1.0e6;
  mout << "# N = " << pn << " ";
  mout << sec << " [SEC] ";
  mout << mups << " [MUPS]" << std::endl;
}
//----------------------------------------------------------------------
