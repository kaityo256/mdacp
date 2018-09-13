//----------------------------------------------------------------------
#include "collision.h"
#include "mpistream.h"
#include "confmaker.h"
//----------------------------------------------------------------------
Collision coll;
//----------------------------------------------------------------------
class BallAdder : public Executor {
private:
  double c[D];
  double r;
  double xv;
  double density;
public:
  BallAdder(double _c[D], double _r, double _xv, double _density) {
    c[X] = _c[X];
    c[Y] = _c[Y];
    c[Z] = _c[Z];
    r = _r;
    xv = _xv;
    density = _density;
  };
  void Add(double x[D], MDUnit *mdu) {
    const double dx = x[X] - c[X];
    const double dy = x[Y] - c[Y];
    const double dz = x[Z] - c[Z];
    const double r2 = dx * dx + dy * dy + dz * dz;
    double v[D] = {xv, 0, 0};
    if (r2  > r * r) return;
    mdu->AddParticle(x, v);
  }
  void Execute(MDUnit *mdu) {
    const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
    const int is = static_cast<int>(2 * r / s);
    const double hs = s * 0.5;
    for (int iz = 0; iz < is; iz++) {
      for (int iy = 0; iy < is; iy++) {
        for (int ix = 0; ix < is; ix++) {
          double x1[D] = {ix * s - r + c[X], iy * s - r + c[Y], iz * s - r + c[Z]};
          double x2[D] = {x1[X] + hs, x1[Y], x1[Z]};
          double x3[D] = {x1[X], x1[Y] + hs, x1[Z]};
          double x4[D] = {x1[X], x1[Y], x1[Z] + hs};
          Add(x1, mdu);
          Add(x2, mdu);
          Add(x3, mdu);
          Add(x4, mdu);
        }
      }
    }
  }
};
//----------------------------------------------------------------------
void
Collision::Run(MDManager *mdm) {
  Parameter *param = mdm->GetParameter();
  const double *L = mdm->GetSystemSize();
  double c1[D] = {L[X] * 0.25, L[Y] * 0.5, L[Z] * 0.5};
  BallAdder ba1(c1, L[Y] * 0.25, 1.0, 0.374);
  double c2[D] = {L[X] * 0.75, L[Y] * 0.5, L[Z] * 0.5};
  BallAdder ba2(c2, L[Y] * 0.25, -1.0, 0.374);
  mdm->ExecuteAll(&ba1);
  mdm->ExecuteAll(&ba2);
  mdm->SaveAsCdviewSequential();
  const int OBSERVE_LOOP = param->GetIntegerDef("ObserveLoop", 100);
  const int TOTAL_LOOP = param->GetIntegerDef("TotalLoop", 1000);
  const bool save_cdview_file = param->GetBooleanDef("SaveCdviewFile", false);
  mdm->ShowSystemInformation();
  for (int i = 0; i < TOTAL_LOOP; i++) {
    mdm->Calculate();
    if (i % OBSERVE_LOOP == 0) {
      mout << mdm->GetSimulationTime();
      mout << " " << mdm->Temperature();
      mout << "# observe" << std::endl;
      if (save_cdview_file) {
        mdm->SaveAsCdviewSequential();
      }
    }
  }
}
//----------------------------------------------------------------------
