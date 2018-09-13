#include <vector>
#include <algorithm>
#include <assert.h>
#include <random>
#include "confmaker.h"
#include "mdunit.h"
//----------------------------------------------------------------------
void
ConfigurationMaker::Execute(MDUnit *mdu) {
  double *L = mdu->GetSystemSize();
  const double density = param->GetDoubleDef("Density", 0.5);
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
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + e;
        x[Y] = static_cast<double>(iy) * s + hs + e;
        x[Z] = static_cast<double>(iz) * s + hs + e;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs + e;
        x[Y] = static_cast<double>(iy) * s + e;
        x[Z] = static_cast<double>(iz) * s + hs + e;
        mdu->AddParticle(x);

        x[X] = static_cast<double>(ix) * s + hs + e;
        x[Y] = static_cast<double>(iy) * s + hs + e;
        x[Z] = static_cast<double>(iz) * s + e;
        mdu->AddParticle(x);
      }
    }
  }
}
//----------------------------------------------------------------------
ExactConfigurationMaker::ExactConfigurationMaker(Parameter *param_, const double *L_) : L(L_) {
  param = param_;
  const double density = param->GetDoubleDef("Density", 0.7);;
  const double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  six = static_cast<int>(L[X] / s) + 1;
  siy = static_cast<int>(L[Y] / s) + 1;
  siz = static_cast<int>(L[Z] / s) + 1;
  unsigned long int n_full = static_cast<unsigned long int>(six * siy * siz * 4);
  unsigned long int n_exact = static_cast<int>(L[X] * L[Y] * L[Z] * density);
  unsigned long int n_delete = n_full - n_exact;
  std::mt19937 mt(1);
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  while (vec_delete.size() < n_delete) {
    unsigned long int v = static_cast<unsigned long int>(n_full * ud(mt));
    if (std::find(vec_delete.begin(), vec_delete.end(), v) == vec_delete.end()) {
      vec_delete.push_back(v);
    }
  }
  std::sort(vec_delete.begin(), vec_delete.end());
}
//----------------------------------------------------------------------
void
ExactConfigurationMaker::Execute(MDUnit *mdu) {
  const double density = param->GetDoubleDef("Density", 0.7);;
  double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  const double hs = s * 0.5;
  const double sx = L[X] / static_cast<double>(six);
  const double sy = L[Y] / static_cast<double>(siy);
  const double sz = L[Z] / static_cast<double>(siz);
  double x[D];
  unsigned long int n = 0;
  unsigned long int n_pos = 0;
  const double e = 0.00000001;
  for (int iz = 0; iz < siz; iz++) {
    for (int iy = 0; iy < siy; iy++) {
      for (int ix = 0; ix < six; ix++) {
        x[X] = static_cast<double>(ix) * sx + e;
        x[Y] = static_cast<double>(iy) * sy + e;
        x[Z] = static_cast<double>(iz) * sz + e;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;

        x[X] = static_cast<double>(ix) * sx + e;
        x[Y] = static_cast<double>(iy) * sy + hs + e;
        x[Z] = static_cast<double>(iz) * sz + hs + e;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;

        x[X] = static_cast<double>(ix) * sx + hs + e;
        x[Y] = static_cast<double>(iy) * sy + e;
        x[Z] = static_cast<double>(iz) * sz + hs + e;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;

        x[X] = static_cast<double>(ix) * sx + hs + e;
        x[Y] = static_cast<double>(iy) * sy + hs + e;
        x[Z] = static_cast<double>(iz) * sz + e;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;
      }
    }
  }
}
//----------------------------------------------------------------------
void
SimpleConfigurationMaker::Execute(MDUnit *mdu) {
  MDRect *myrect = mdu->GetRect();
  double w[D];
  w[X] = myrect->GetWidth(X);
  w[Y] = myrect->GetWidth(Y);
  w[Z] = myrect->GetWidth(Z);
  const double e = 0.0000001;
  assert((w[X] - w[Y]) / w[X] < e && (w[X] - w[Y]) / w[X] > -e);
  assert((w[X] - w[Z]) / w[X] < e && (w[X] - w[Z]) / w[X] > -e);
  const double density = param->GetDoubleDef("Density", 0.5);
  double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  int sx = static_cast<unsigned long int>(w[X] / s);
  int sy = static_cast<unsigned long int>(w[Y] / s);
  int sz = static_cast<unsigned long int>(w[Z] / s);

  s = w[X] / sx;
  double hs = s * 0.5;

  double x[D];
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        x[X] = myrect->s[X] + static_cast<double>(ix) * s + e;
        x[Y] = myrect->s[Y] + static_cast<double>(iy) * s + e;
        x[Z] = myrect->s[Z] + static_cast<double>(iz) * s + e;
        mdu->AddParticle(x);

        x[X] = myrect->s[X] + static_cast<double>(ix) * s + e;
        x[Y] = myrect->s[Y] + static_cast<double>(iy) * s + hs + e;
        x[Z] = myrect->s[Z] + static_cast<double>(iz) * s + hs + e;
        mdu->AddParticle(x);

        x[X] = myrect->s[X] + static_cast<double>(ix) * s + hs + e;
        x[Y] = myrect->s[Y] + static_cast<double>(iy) * s + e;
        x[Z] = myrect->s[Z] + static_cast<double>(iz) * s + hs + e;
        mdu->AddParticle(x);

        x[X] = myrect->s[X] + static_cast<double>(ix) * s + hs + e;
        x[Y] = myrect->s[Y] + static_cast<double>(iy) * s + hs + e;
        x[Z] = myrect->s[Z] + static_cast<double>(iz) * s + e;
        mdu->AddParticle(x);
      }
    }
  }
}
//----------------------------------------------------------------------
void
ExtractConfigurationMaker::Execute(MDUnit *mdu) {
  MDRect *myrect = mdu->GetRect();
  double w[D];
  w[X] = myrect->GetWidth(X);
  w[Y] = myrect->GetWidth(Y);
  w[Z] = myrect->GetWidth(Z);
  const double e = 0.00000001;
  assert((w[X] - w[Y]) / w[X] < e && (w[X] - w[Y]) / w[X] > -e);
  assert((w[X] - w[Z]) / w[X] < e && (w[X] - w[Z]) / w[X] > -e);
  double s = 1.0 / pow(density * 0.25, 1.0 / 3.0);
  int sx = static_cast<unsigned long int>(w[X] / s) + 1;
  int sy = static_cast<unsigned long int>(w[Y] / s) + 1;
  int sz = static_cast<unsigned long int>(w[Z] / s) + 1;
  const int n_full = sx * sy * sz * 4;
  const int n_exact = static_cast<int>(w[X] * w[Y] * w[Z] * density);
  const int n_delete = n_full - n_exact;
  std::vector<int> vec_delete;
  std::mt19937 mt(mdu->GetID());
  std::uniform_real_distribution<double> ud(0.0, 1.0);
  while ((int)vec_delete.size() < n_delete) {
    int v = static_cast<unsigned long int>(n_full * ud(mt));
    if (std::find(vec_delete.begin(), vec_delete.end(), v) == vec_delete.end()) {
      vec_delete.push_back(v);
    }
  }
  std::sort(vec_delete.begin(), vec_delete.end());

  s = w[X] / sx;
  double hs = s * 0.5;

  double x[D];
  int n = 0;
  int n_pos = 0;
  for (int iz = 0; iz < sz; iz++) {
    for (int iy = 0; iy < sy; iy++) {
      for (int ix = 0; ix < sx; ix++) {
        x[X] = myrect->s[X] + static_cast<double>(ix) * s;
        x[Y] = myrect->s[Y] + static_cast<double>(iy) * s;
        x[Z] = myrect->s[Z] + static_cast<double>(iz) * s;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;

        x[X] = myrect->s[X] + static_cast<double>(ix) * s;
        x[Y] = myrect->s[Y] + static_cast<double>(iy) * s + hs;
        x[Z] = myrect->s[Z] + static_cast<double>(iz) * s + hs;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;

        x[X] = myrect->s[X] + static_cast<double>(ix) * s + hs;
        x[Y] = myrect->s[Y] + static_cast<double>(iy) * s;
        x[Z] = myrect->s[Z] + static_cast<double>(iz) * s + hs;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;

        x[X] = myrect->s[X] + static_cast<double>(ix) * s + hs;
        x[Y] = myrect->s[Y] + static_cast<double>(iy) * s + hs;
        x[Z] = myrect->s[Z] + static_cast<double>(iz) * s;
        if (vec_delete[n_pos] == n) {
          n_pos++;
        } else {
          mdu->AddParticle(x);
        }
        n++;
      }
    }
  }
}
//----------------------------------------------------------------------
