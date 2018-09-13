//----------------------------------------------------------------------
#ifndef mdrect_h
#define mdrect_h
//----------------------------------------------------------------------
#include <stdio.h>
#include "mdconfig.h"
#include "simulationinfo.h"
//----------------------------------------------------------------------
class MDRect {
public:
  double s[D], e[D];
  MDRect(void) {
    for (int i = 0; i < D; i++) {
      s[i] = 0.0;
      e[i] = 0.0;
    }
  };
  MDRect(double s2[D], double e2[D]) {
    for (int i = 0; i < D; i++) {
      s[i] = s2[i];
      e[i] = e2[i];
    }
  };
  MDRect(double sx, double ex, double sy, double ey, double sz, double ez) {
    s[X] = sx;
    e[X] = ex;
    s[Y] = sy;
    e[Y] = ey;
    s[Z] = sz;
    e[Z] = ez;
  };
  double GetWidth(int dir) {return e[dir] - s[dir];}

  bool IsInsideEdge(int dir, double q[D], SimulationInfo *sinfo);
  bool IsInside(double q[D]) {
    return (s[X] <= q[X] && q[X] < e[X] && s[Y] <= q[Y] && q[Y] < e[Y] && s[Z] <= q[Z] &&
            q[Z] < e[Z]);
  };
  void Show(void) {
    printf("(%f,%f,%f)->(%f,%f,%f)\n", s[X], s[Y], s[Z], e[X], e[Y], e[Z]);
  };
  double * GetStartPosition(void) {return s;}
  bool IsOverBoundary(int dir, double q[D]);
  void ChangeScale(double alpha);
};
//----------------------------------------------------------------------
std::ostream& operator<<(std::ostream& os, const MDRect& rect);
//----------------------------------------------------------------------
#endif
