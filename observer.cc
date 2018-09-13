//----------------------------------------------------------------------
#include "observer.h"
//----------------------------------------------------------------------
double
KineticEnergyObserver::Observe(Variables *vars, MeshList *) {
  double e = 0;
  const int pn = vars->GetParticleNumber();
  double (*p)[D] = vars->p;
  for (int i = 0; i < pn; i++) {
    for (int d = 0; d < D; d++) {
      e += 0.5 * p[i][d] * p[i][d];
    }
  }
  return e;
}
//----------------------------------------------------------------------
double
PotentialEnergyObserver::Observe(Variables *vars, MeshList *mesh) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  double (*q)[D] = vars->q;
  const double C2 = vars->GetC2();
  const double C0 = vars->GetC0();
  const int pn = vars->GetParticleNumber();

  double energy = 0.0;

  const int s = mesh->GetPairNumber();
  int *key_particles = mesh->GetKeyParticles();
  int *partner_particles = mesh->GetPartnerParticles();
  for (int k = 0; k < s; k++) {
    int i = key_particles[k];
    int j = partner_particles[k];
    double dx = q[i][X] - q[j][X];
    double dy = q[i][Y] - q[j][Y];
    double dz = q[i][Z] - q[j][Z];
    const double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2) continue;
    double e = 4.0 * (1.0 / (r2 * r2 * r2 * r2 * r2 * r2) - 1.0 / (r2 * r2 * r2) + C2 * r2 + C0);
    if (i >= pn || j >= pn ) {
      e *= 0.5;
    }
    energy += e;
  }
#ifdef FX10
  energy *= 0.5;
#endif
  return energy;
}
//----------------------------------------------------------------------
double
VirialObserver::Observe(Variables *vars, MeshList *mesh) {
  double phi = 0.0;
  double (*q)[D] = vars->q;
  const int pn = vars->GetParticleNumber();
  const double C_2 = vars->GetC2();
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const int s = mesh->GetPairNumber();
  int *key_particles = mesh->GetKeyParticles();
  int *partner_particles = mesh->GetPartnerParticles();
  for (int k = 0; k < s; k++) {
    int i = key_particles[k];
    int j = partner_particles[k];
    double dx = q[j][X] - q[i][X];
    double dy = q[j][Y] - q[i][Y];
    double dz = q[j][Z] - q[i][Z];
    const double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2) continue;
    const double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C_2 * 8.0) * r2;
    if ( i >= pn || j >= pn) {
      df *= 0.5;
    }
    phi += df;
  }
#ifdef FX10
  phi *= 0.5;
#endif
  return phi / 3.0;
}
//----------------------------------------------------------------------
