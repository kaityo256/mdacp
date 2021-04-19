//---------------------------------------------------------------------
#ifndef mdconfig_h
#define mdconfig_h
//---------------------------------------------------------------------
#include <iostream>
#include <stdlib.h>
#include <math.h>

const int D = 4;
const int X = 0, Y = 1, Z = 2;
const int N = 1000000;
//const int PAIRLIST_SIZE = N * 80;
const int PAIRLIST_SIZE = N*30;

const double CUTOFF_LENGTH = 3.0;
//const double CUTOFF_LENGTH = 2.5;

//---------------------------------------------------------------------------
extern const char *MDACP_VERSION;
//---------------------------------------------------------------------------
#define show_error(MSG) { mout << "# Error at " << __FILE__ <<":" << __LINE__ << std::endl;mout << MSG << std::endl;}
#define show_warning(MSG) { mout << "# Warning at " << __FILE__ <<":" << __LINE__ << std::endl;mout << MSG << std::endl;}
//---------------------------------------------------------------------------
const int MAX_DIR = 6;
enum DIRECTION {D_LEFT, D_RIGHT, D_BACK, D_FORWARD, D_DOWN, D_UP};
const int OppositeDir[MAX_DIR] = {D_RIGHT, D_LEFT, D_FORWARD, D_BACK, D_UP, D_DOWN};
enum HEATBATH_TYPE {HT_NOSEHOOVER, HT_LANGEVIN};
//---------------------------------------------------------------------------
class Direction {
private:
  static const char *name_str[MAX_DIR];
public:
  static const char *Name(int dir) {
    return name_str[dir];
  }
};
//---------------------------------------------------------------------------
struct ParticleInfo {
  double q[D];
  double p[D];
  int type;
  ParticleInfo() {
    q[X] = 0.0;
    q[Y] = 0.0;
    q[Z] = 0.0;
    p[X] = 0.0;
    p[Y] = 0.0;
    p[Z] = 0.0;
    type = 0;
  };
  ParticleInfo(double qx, double qy, double qz, double px, double py, double pz, int t) {
    q[X] = qx;
    q[Y] = qy;
    q[Z] = qz;
    p[X] = px;
    p[Y] = py;
    p[Z] = pz;
    type = t;
  }
};
//---------------------------------------------------------------------------
#endif
//---------------------------------------------------------------------
