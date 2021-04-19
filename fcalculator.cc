//----------------------------------------------------------------------
#include <iostream>
#include <fstream>
#include <random>
#include "fcalculator.h"
//----------------------------------------------------------------------
#ifdef FX10
#include "fj_tool/fapp.h"
//#include "simd_fx10.h"
#endif
#ifdef AVX2
#include <immintrin.h>
typedef double v4df __attribute__((vector_size(32)));
void
put(v4df r) {
  double *a = (double*)(&r);
  printf("%.10f %.10f %.10f %.10f\n", a[0], a[1], a[2], a[3]);
}
#endif
//----------------------------------------------------------------------
void
ForceCalculator::UpdatePositionHalf(Variables *vars, SimulationInfo *sinfo) {
  const double dt2 = sinfo->TimeStep * 0.5;
  const int pn = vars->GetParticleNumber();
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  for (int i = 0; i < pn; i++) {
    q[i][X] += p[i][X] * dt2;
    q[i][Y] += p[i][Y] * dt2;
    q[i][Z] += p[i][Z] * dt2;
  }
}
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForce(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {

#ifdef FX10
fapp_start("force", 0,0);
  CalculateForceReactless(vars,mesh,sinfo);
fapp_stop("force", 0,0);
//  CalculateForceReactlessSIMD(vars,mesh,sinfo);
  //CalculateForceReactlessSIMD_errsafe(vars, mesh, sinfo);
#elif AVX2
  CalculateForceAVX2(vars, mesh, sinfo);
#else
  CalculateForceNext(vars, mesh, sinfo);
  //CalculateForceReactless(vars,mesh,sinfo);
  //CalculateForceBruteforce(vars,sinfo);
  //CalculateForceSorted(vars,mesh,sinfo);
  //CalculateForcePair(vars,mesh,sinfo);
  //CalculateForceUnroll(vars,mesh,sinfo);
#endif
}
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceBruteforce(Variables *vars, SimulationInfo *sinfo) {
  const double dt = sinfo->TimeStep;
  const int pn = vars->GetParticleNumber();
  const int tn = vars->GetTotalParticleNumber();
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  for (int i = 0; i < pn; i++) {
    for (int j = i + 1; j < tn; j++) {
      const double dx = q[j][X] - q[i][X];
      const double dy = q[j][Y] - q[i][Y];
      const double dz = q[j][Z] - q[i][Z];
      const double r2 = (dx * dx + dy * dy + dz * dz);
      if (r2 > CL2) continue;
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      p[i][X] += df * dx;
      p[i][Y] += df * dy;
      p[i][Z] += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
  }
}
//----------------------------------------------------------------------
// CalculateForce  (Optimized for IBM POWER )
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceUnroll(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double C2_8 = C2 * 8.0;
  const double dt = sinfo->TimeStep;

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int pn = vars->GetParticleNumber();

  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = mesh->GetPartnerNumber(i);
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    for (int k = kp; k < (kp + np - 1); k += 2) {
      const int j_a = sorted_list[k];
      const int j_b = sorted_list[k + 1];
      const double dx_a = q[j_a][X] - qx_key;
      const double dy_a = q[j_a][Y] - qy_key;
      const double dz_a = q[j_a][Z] - qz_key;

      const double dx_b = q[j_b][X] - qx_key;
      const double dy_b = q[j_b][Y] - qy_key;
      const double dz_b = q[j_b][Z] - qz_key;

      const double r2_a = (dx_a * dx_a + dy_a * dy_a + dz_a * dz_a);
      const double r2_b = (dx_b * dx_b + dy_b * dy_b + dz_b * dz_b);

      const double r6_a = r2_a * r2_a * r2_a;
      const double r6_b = r2_b * r2_b * r2_b;

      const double r14_a = r6_a * r6_a * r2_a;
      const double r14_b = r6_b * r6_b * r2_b;

      const double r14_inv_ab = 1.0 / (r14_a * r14_b);
      double df_a = (24.0 * r6_a - 48.0) * r14_inv_ab * r14_b;
      double df_b = (24.0 * r6_b - 48.0) * r14_inv_ab * r14_a;

      df_a = (df_a + C2_8) * dt;
      df_b = (df_b + C2_8) * dt;

      if (r2_a > CL2) {
        df_a = 0.0;
      }

      if (r2_b > CL2) {
        df_b = 0.0;
      }

      pfx += df_a * dx_a;
      pfy += df_a * dy_a;
      pfz += df_a * dz_a;

      pfx += df_b * dx_b;
      pfy += df_b * dy_b;
      pfz += df_b * dz_b;

      p[j_a][X] -= df_a * dx_a;
      p[j_a][Y] -= df_a * dy_a;
      p[j_a][Z] -= df_a * dz_a;

      p[j_b][X] -= df_b * dx_b;
      p[j_b][Y] -= df_b * dy_b;
      p[j_b][Z] -= df_b * dz_b;
    }
    if (np % 2 == 1) {
      int j = sorted_list[kp + np - 1];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      if (r2 > CL2) {
        df = 0.0;
      }
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//----------------------------------------------------------------------
// Calculate Force using sorted list.
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceSorted(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double dt = sinfo->TimeStep;
  const int pn = vars->GetParticleNumber();

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = mesh->GetPartnerNumber(i);
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      if (r2 > CL2) {
        df = 0.0;
      }
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//----------------------------------------------------------------------
// Calculate Force (Optimized for Intel)
// Calculate Next Pair on previous loop
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceNext(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {

  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2() * 8.0;
  const double dt = sinfo->TimeStep;

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int pn = vars->GetParticleNumber();

  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    int ja = sorted_list[kp];
    double dxa = q[ja][X] - qx_key;
    double dya = q[ja][Y] - qy_key;
    double dza = q[ja][Z] - qz_key;
    double df = 0.0;
    double dxb = 0.0, dyb = 0.0, dzb = 0.0;
    int jb = 0;

    const int np = mesh->GetPartnerNumber(i);
    for (int k = kp; k < np + kp; k++) {

      const double dx = dxa;
      const double dy = dya;
      const double dz = dza;
      double r2 = (dx * dx + dy * dy + dz * dz);
      const int j = ja;
      ja = sorted_list[k + 1];
      dxa = q[ja][X] - qx_key;
      dya = q[ja][Y] - qy_key;
      dza = q[ja][Z] - qz_key;
      if (r2 > CL2)continue;
      pfx += df * dxb;
      pfy += df * dyb;
      pfz += df * dzb;
      p[jb][X] -= df * dxb;
      p[jb][Y] -= df * dyb;
      p[jb][Z] -= df * dzb;
      const double r6 = r2 * r2 * r2;
      df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2) * dt;
      jb = j;
      dxb = dx;
      dyb = dy;
      dzb = dz;
    }
    p[jb][X] -= df * dxb;
    p[jb][Y] -= df * dyb;
    p[jb][Z] -= df * dzb;
    p[i][X] += pfx + df * dxb;
    p[i][Y] += pfy + df * dyb;
    p[i][Z] += pfz + df * dzb;
  }
}
//----------------------------------------------------------------------
#ifdef AVX2
void
ForceCalculator::CalculateForceAVX2(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {

  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2() * 8.0;
  const double dt = sinfo->TimeStep;
  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int pn = vars->GetParticleNumber();
  const int *sorted_list = mesh->GetSortedList();

  const v4df vzero = _mm256_set_pd(0, 0, 0, 0);
  const v4df vcl2 = _mm256_set_pd(CL2, CL2, CL2, CL2);
  const v4df vc24 = _mm256_set_pd(24 * dt, 24 * dt, 24 * dt, 24 * dt);
  const v4df vc48 = _mm256_set_pd(48 * dt, 48 * dt, 48 * dt, 48 * dt);
  const v4df vc2 = _mm256_set1_pd(C2 * dt);

  for (int i = 0; i < pn; i++) {
    const int np = mesh->GetPartnerNumber(i);
    const v4df vqi = _mm256_load_pd((double*)(q + i));
    v4df vpf = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    const int kp = mesh->GetKeyPointer(i);
    int ja_1 = sorted_list[kp];
    int ja_2 = sorted_list[kp + 1];
    int ja_3 = sorted_list[kp + 2];
    int ja_4 = sorted_list[kp + 3];
    v4df vqj_1 = _mm256_load_pd((double*)(q + ja_1));
    v4df vdqa_1 = vqj_1 - vqi;
    v4df vqj_2 = _mm256_load_pd((double*)(q + ja_2));
    v4df vdqa_2 = vqj_2 - vqi;
    v4df vqj_3 = _mm256_load_pd((double*)(q + ja_3));
    v4df vdqa_3 = vqj_3 - vqi;
    v4df vqj_4 = _mm256_load_pd((double*)(q + ja_4));
    v4df vdqa_4 = vqj_4 - vqi;

    v4df vdf = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);

    v4df vdqb_1 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    v4df vdqb_2 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    v4df vdqb_3 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);
    v4df vdqb_4 = _mm256_set_pd(0.0, 0.0, 0.0, 0.0);

    int jb_1 = 0, jb_2 = 0, jb_3 = 0, jb_4 = 0;
    for (int k = 0; k < (np / 4) * 4; k += 4) {
      const int j_1 = ja_1;
      const int j_2 = ja_2;
      const int j_3 = ja_3;
      const int j_4 = ja_4;
      v4df vdq_1 = vdqa_1;
      v4df vdq_2 = vdqa_2;
      v4df vdq_3 = vdqa_3;
      v4df vdq_4 = vdqa_4;

      ja_1 = sorted_list[kp + k + 4];
      ja_2 = sorted_list[kp + k + 5];
      ja_3 = sorted_list[kp + k + 6];
      ja_4 = sorted_list[kp + k + 7];

      v4df tmp0 = _mm256_unpacklo_pd(vdq_1, vdq_2);
      v4df tmp1 = _mm256_unpackhi_pd(vdq_1, vdq_2);
      v4df tmp2 = _mm256_unpacklo_pd(vdq_3, vdq_4);
      v4df tmp3 = _mm256_unpackhi_pd(vdq_3, vdq_4);

      v4df vdx = _mm256_permute2f128_pd(tmp0, tmp2, 0x20);
      v4df vdy = _mm256_permute2f128_pd(tmp1, tmp3, 0x20);
      v4df vdz = _mm256_permute2f128_pd(tmp0, tmp2, 0x31);

      v4df vdf_1 = _mm256_permute4x64_pd(vdf, 0);
      v4df vdf_2 = _mm256_permute4x64_pd(vdf, 85);
      v4df vdf_3 = _mm256_permute4x64_pd(vdf, 170);
      v4df vdf_4 = _mm256_permute4x64_pd(vdf, 255);

      vqj_1 = _mm256_load_pd((double*)(q + ja_1));
      vdqa_1 = vqj_1 - vqi;
      vpf += vdf_1 * vdqb_1;

      v4df vpjb_1 = _mm256_load_pd((double*)(p + jb_1));
      vpjb_1 -= vdf_1 * vdqb_1;
      _mm256_store_pd((double*)(p + jb_1), vpjb_1);

      vqj_2 = _mm256_load_pd((double*)(q + ja_2));
      vdqa_2 = vqj_2 - vqi;
      vpf += vdf_2 * vdqb_2;

      v4df vpjb_2 = _mm256_load_pd((double*)(p + jb_2));
      vpjb_2 -= vdf_2 * vdqb_2;
      _mm256_store_pd((double*)(p + jb_2), vpjb_2);

      vqj_3 = _mm256_load_pd((double*)(q + ja_3));
      vdqa_3 = vqj_3 - vqi;
      vpf += vdf_3 * vdqb_3;

      v4df vpjb_3 = _mm256_load_pd((double*)(p + jb_3));
      vpjb_3 -= vdf_3 * vdqb_3;
      _mm256_store_pd((double*)(p + jb_3), vpjb_3);

      vqj_4 = _mm256_load_pd((double*)(q + ja_4));
      vdqa_4 = vqj_4 - vqi;
      vpf += vdf_4 * vdqb_4;

      v4df vpjb_4 = _mm256_load_pd((double*)(p + jb_4));
      vpjb_4 -= vdf_4 * vdqb_4;
      _mm256_store_pd((double*)(p + jb_4), vpjb_4);

      v4df vr2 = vdx * vdx + vdy * vdy + vdz * vdz;
      v4df vr6 = vr2 * vr2 * vr2;
      vdf = (vc24 * vr6 - vc48) / (vr6 * vr6 * vr2) + vc2;
      v4df mask = vcl2 - vr2;
      vdf = _mm256_blendv_pd(vdf, vzero, mask);

      jb_1 = j_1;
      jb_2 = j_2;
      jb_3 = j_3;
      jb_4 = j_4;
      vdqb_1 = vdq_1;
      vdqb_2 = vdq_2;
      vdqb_3 = vdq_3;
      vdqb_4 = vdq_4;
    }
    v4df vdf_1 = _mm256_permute4x64_pd(vdf, 0);
    v4df vdf_2 = _mm256_permute4x64_pd(vdf, 85);
    v4df vdf_3 = _mm256_permute4x64_pd(vdf, 170);
    v4df vdf_4 = _mm256_permute4x64_pd(vdf, 255);

    v4df vpjb_1 = _mm256_load_pd((double*)(p + jb_1));
    vpjb_1 -= vdf_1 * vdqb_1;
    _mm256_store_pd((double*)(p + jb_1), vpjb_1);

    v4df vpjb_2 = _mm256_load_pd((double*)(p + jb_2));
    vpjb_2 -= vdf_2 * vdqb_2;
    _mm256_store_pd((double*)(p + jb_2), vpjb_2);

    v4df vpjb_3 = _mm256_load_pd((double*)(p + jb_3));
    vpjb_3 -= vdf_3 * vdqb_3;
    _mm256_store_pd((double*)(p + jb_3), vpjb_3);

    v4df vpjb_4 = _mm256_load_pd((double*)(p + jb_4));
    vpjb_4 -= vdf_4 * vdqb_4;
    _mm256_store_pd((double*)(p + jb_4), vpjb_4);

    v4df vpi = _mm256_load_pd((double*)(p + i));
    vpf += vdf_1 * vdqb_1;
    vpf += vdf_2 * vdqb_2;
    vpf += vdf_3 * vdqb_3;
    vpf += vdf_4 * vdqb_4;
    vpi += vpf;
    _mm256_store_pd((double*)(p + i), vpi);
    const double qix = q[i][X];
    const double qiy = q[i][Y];
    const double qiz = q[i][Z];
    double pfx = 0.0;
    double pfy = 0.0;
    double pfz = 0.0;
    for (int k = (np / 4) * 4; k < np; k++) {
      const int j = sorted_list[k + kp];
      double dx = q[j][X] - qix;
      double dy = q[j][Y] - qiy;
      double dz = q[j][Z] - qiz;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2) * dt;
      if (r2 > CL2) df = 0.0;
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
      p[j][X] -= df * dx;
      p[j][Y] -= df * dy;
      p[j][Z] -= df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
#endif
//----------------------------------------------------------------------
// Calculate Force without optimization
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForcePair(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double dt = sinfo->TimeStep;

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;

  int *key_index = mesh->GetKeyParticles();
  int *partner_index = mesh->GetPartnerParticles();

  const int number_of_pairs = mesh->GetPairNumber();
  for (int k = 0; k < number_of_pairs; k++) {
    int i = key_index[k];
    int j = partner_index[k];
    double dx = q[j][X] - q[i][X];
    double dy = q[j][Y] - q[i][Y];
    double dz = q[j][Z] - q[i][Z];
    double r2 = (dx * dx + dy * dy + dz * dz);
    if (r2 > CL2)continue;
    double r6 = r2 * r2 * r2;
    double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
    p[i][X] += df * dx;
    p[i][Y] += df * dy;
    p[i][Z] += df * dz;
    p[j][X] -= df * dx;
    p[j][Y] -= df * dy;
    p[j][Z] -= df * dz;
  }
}
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceReactless(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double CL2 = CUTOFF_LENGTH * CUTOFF_LENGTH;
  const double C2 = vars->GetC2();
  const double dt = sinfo->TimeStep;
  const int pn = vars->GetTotalParticleNumber();

  double (*q)[D] = vars->q;
  double (*p)[D] = vars->p;
  const int *sorted_list = mesh->GetSortedList();

  for (int i = 0; i < pn; i++) {
    const double qx_key = q[i][X];
    const double qy_key = q[i][Y];
    const double qz_key = q[i][Z];
    const int np = mesh->GetPartnerNumber(i);
    double pfx = 0;
    double pfy = 0;
    double pfz = 0;
    const int kp = mesh->GetKeyPointer(i);
    for (int k = 0; k < np; k++) {
      const int j = sorted_list[kp + k];
      double dx = q[j][X] - qx_key;
      double dy = q[j][Y] - qy_key;
      double dz = q[j][Z] - qz_key;
      double r2 = (dx * dx + dy * dy + dz * dz);
      double r6 = r2 * r2 * r2;
      double df = ((24.0 * r6 - 48.0) / (r6 * r6 * r2) + C2 * 8.0) * dt;
      if (r2 > CL2) {
        df = 0.0;
      }
      pfx += df * dx;
      pfy += df * dy;
      pfz += df * dz;
    }
    p[i][X] += pfx;
    p[i][Y] += pfy;
    p[i][Z] += pfz;
  }
}
//----------------------------------------------------------------------
#if 0
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceReactlessSIMD(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double C2 = vars->GetC2();
  const double dt = sinfo->TimeStep;

  double (*p)[D] = vars->p;
  const int pn = vars->GetTotalParticleNumber();
  static double q[N][D];
  memcpy((void*)q, (void*)vars->q, sizeof(double) * 3 * pn);

  const int *sorted_list = mesh->GetSortedList();
  int *key_pointer = mesh->GetKeyPointerP();
  int *number_of_partners = mesh->GetNumberOfPartners();
  const double Rcf = CUTOFF_LENGTH;
  const double Rcf2 = Rcf * Rcf;
  const double c1 = dt * (48 * pow(Rcf, -14) - 24 * pow(Rcf, -8));
  v2df c1_v2(c1, c1);
  v2df _Rcf2v2(1.0 / Rcf2, 1.0 / Rcf2);


  for (int i = 0; i < pn; i++) {
    v2df qxi(q[i][X], q[i][X]);
    v2df qyi(q[i][Y], q[i][Y]);
    v2df qzi(q[i][Z], q[i][Z]);
    const int kp = key_pointer[i];
    const int np = number_of_partners[i];
    v2df vxiv2(0.0, 0.0), vyiv2(0.0, 0.0), vziv2(0.0, 0.0);
    int lj = kp;
    for (; lj < kp + np - 3; lj += 4) {
      int ja = sorted_list[lj];
      int jb = sorted_list[lj + 1];
      int jc = sorted_list[lj + 2];
      int jd = sorted_list[lj + 3];
      v2df dxa(q[ja][X], q[jb][X]);
      v2df dya(q[ja][Y], q[jb][Y]);
      v2df dza(q[ja][Z], q[jb][Z]);
      v2df dxc(q[jc][X], q[jd][X]);
      v2df dyc(q[jc][Y], q[jd][Y]);
      v2df dzc(q[jc][Z], q[jd][Z]);
      dxa -= qxi;
      dxc -= qxi;
      dya -= qyi;
      dyc -= qyi;
      dza -= qzi;
      dzc -= qzi;
      v2df RRa = dxa * dxa + dya * dya + dza * dza;
      v2df RRc = dxc * dxc + dyc * dyc + dzc * dzc;
      v2df mask_a = v2df(1.0, 1.0) - (RRa * _Rcf2v2).floor();
      v2df mask_c = v2df(1.0, 1.0) - (RRc * _Rcf2v2).floor();
      v2df R4a = RRa * RRa;
      v2df R4c = RRc * RRc;
      v2df R6a = R4a * RRa;
      v2df R6c = R4c * RRc;
      v2df R8a = R4a * R4a;
      v2df R8c = R4c * R4c;
      v2df dv_a = (v2df(48 * dt, 48 * dt) - v2df(24 * dt, 24 * dt) * R6a) / (R8a * R6a) - c1_v2;
      v2df dv_c = (v2df(48 * dt, 48 * dt) - v2df(24 * dt, 24 * dt) * R6c) / (R8c * R6c) - c1_v2;
      dxa = mask_a * dxa;
      dya = mask_a * dya;
      dza = mask_a * dza;
      dxc = mask_c * dxc;
      dyc = mask_c * dyc;
      dzc = mask_c * dzc;

      vxiv2 += dxa * dv_a + dxc * dv_c;
      vyiv2 += dya * dv_a + dyc * dv_c;
      vziv2 += dza * dv_a + dzc * dv_c;
    }
    p[i][X] -= vxiv2[0] + vxiv2[1];
    p[i][Y] -= vyiv2[0] + vyiv2[1];
    p[i][Z] -= vziv2[0] + vziv2[1];
    for ( ; lj < kp + np; lj++) {
      int j = sorted_list[lj];
      double dx = q[j][X] - q[i][X];
      double dy = q[j][Y] - q[i][Y];
      double dz = q[j][Z] - q[i][Z];
      double R2 = dx * dx + dy * dy + dz * dz;
      if (R2 < Rcf2) {
        double R4 = R2 * R2;
        double R6 = R4 * R2;
        double R8 = R4 * R4;
        double dv = (48 * dt - 24 * dt * R6) / (R8 * R6) - c1;
        p[i][X] -= dx * dv;
        p[i][Y] -= dy * dv;
        p[i][Z] -= dz * dv;
      }
    }
  }

}
//----------------------------------------------------------------------
void
ForceCalculator::CalculateForceReactlessSIMD_errsafe(Variables *vars, MeshList *mesh, SimulationInfo *sinfo) {
  const double C2 = vars->GetC2();
  const double dt = sinfo->TimeStep;

  double (*p)[D] = vars->p;
  const int pn = vars->GetTotalParticleNumber();
  //double (*q)[D] = vars->q;
  thread_local static double q[N][D];
  memcpy((void*)q, (void*)vars->q, sizeof(double) * 4 * pn);

  const int *sorted_list = mesh->GetSortedList();
  int *key_pointer = mesh->GetKeyPointerP();
  int *number_of_partners = mesh->GetNumberOfPartners();

  const double Rcf = CUTOFF_LENGTH;
  const double Rcf2 = Rcf * Rcf;
  const double c1 = dt * (48 * pow(Rcf, -14) - 24 * pow(Rcf, -8));
  v2df c1_v2(c1, c1);
  v2df _Rcf2v2(1.0 / Rcf2, 1.0 / Rcf2);


  for (int i = 0; i < pn; i++) {
    v2df qxi(q[i][X], q[i][X]);
    v2df qyi(q[i][Y], q[i][Y]);
    v2df qzi(q[i][Z], q[i][Z]);
    const int kp = key_pointer[i];
    const int np = number_of_partners[i];
    v2df vxiv2(0.0, 0.0), vyiv2(0.0, 0.0), vziv2(0.0, 0.0);
    int lj = kp;
    for (; lj < kp + np - 3; lj += 4) {
      int ja = sorted_list[lj];
      int jb = sorted_list[lj + 1];
      int jc = sorted_list[lj + 2];
      int jd = sorted_list[lj + 3];
      v2df dxa(q[ja][X], q[jb][X]);
      v2df dya(q[ja][Y], q[jb][Y]);
      v2df dza(q[ja][Z], q[jb][Z]);
      v2df dxc(q[jc][X], q[jd][X]);
      v2df dyc(q[jc][Y], q[jd][Y]);
      v2df dzc(q[jc][Z], q[jd][Z]);
      dxa -= qxi;
      dxc -= qxi;
      dya -= qyi;
      dyc -= qyi;
      dza -= qzi;
      dzc -= qzi;
      v2df RRa = dxa * dxa + dya * dya + dza * dza;
      v2df RRc = dxc * dxc + dyc * dyc + dzc * dzc;
      v2df _R2a = RRa.inv();
      v2df _R2c = RRc.inv();
      v2df mask_a = (RRa * _Rcf2v2).floor();
      v2df mask_c = (RRc * _Rcf2v2).floor();
      v2df Aa = v2df(8 * 48 * dt, 8 * 48 * dt) - v2df(7 * 48 * dt, 7 * 48 * dt) * RRa * _R2a;
      v2df Ac = v2df(8 * 48 * dt, 8 * 48 * dt) - v2df(7 * 48 * dt, 7 * 48 * dt) * RRc * _R2c;
      v2df Ba = v2df(5 * 24 * dt, 5 * 24 * dt) - v2df(4 * 24 * dt, 4 * 24 * dt) * RRa * _R2a;
      v2df Bc = v2df(5 * 24 * dt, 5 * 24 * dt) - v2df(4 * 24 * dt, 4 * 24 * dt) * RRc * _R2c;
      v2df _R4a = _R2a * _R2a;
      v2df _R4c = _R2c * _R2c;
      v2df _R6a = _R4a * _R2a;
      v2df _R6c = _R4c * _R2c;
      v2df _R8a = _R4a * _R4a;
      v2df _R8c = _R4c * _R4c;
      v2df dv_a = Aa * _R8a * _R6a - Ba * _R8a - c1_v2;
      v2df dv_c = Ac * _R8c * _R6c - Bc * _R8c - c1_v2;
      dxa = dxa - dxa * mask_a;
      dya = dya - dya * mask_a;
      dza = dza - dza * mask_a;
      dxc = dxc - dxc * mask_c;
      dyc = dyc - dyc * mask_c;
      dzc = dzc - dzc * mask_c;

      vxiv2 += dxa * dv_a + dxc * dv_c;
      vyiv2 += dya * dv_a + dyc * dv_c;
      vziv2 += dza * dv_a + dzc * dv_c;
    }
    p[i][X] -= vxiv2[0] + vxiv2[1];
    p[i][Y] -= vyiv2[0] + vyiv2[1];
    p[i][Z] -= vziv2[0] + vziv2[1];
    //printf("%f %f %f\n",p[i][X],p[i][Y],p[i][Z]);
    for ( ; lj < kp + np; lj++) {
      int j = sorted_list[lj];
      double dx = q[j][X] - q[i][X];
      double dy = q[j][Y] - q[i][Y];
      double dz = q[j][Z] - q[i][Z];
      double R2 = dx * dx + dy * dy + dz * dz;
      if (R2 < Rcf2) {
        double R4 = R2 * R2;
        double R6 = R4 * R2;
        double R8 = R4 * R4;
        double dv = (48 * dt - 24 * dt * R6) / (R8 * R6) - c1;
        p[i][X] -= dx * dv;
        p[i][Y] -= dy * dv;
        p[i][Z] -= dz * dv;
      }
    }
  }

}
//----------------------------------------------------------------------
#endif
//----------------------------------------------------------------------

//----------------------------------------------------------------------
// For Heatbath
//----------------------------------------------------------------------
void
ForceCalculator::HeatbathZeta(Variables *vars, double current_temperature, SimulationInfo *sinfo) {
  const double dt2 = sinfo->TimeStep * 0.5;
  const double tau = sinfo->HeatbathTau;
  double t1 = (current_temperature - sinfo->AimedTemperature) / (tau * tau);
  vars->Zeta += t1 * dt2;
}
//----------------------------------------------------------------------
void
ForceCalculator::HeatbathMomenta(Variables *vars, SimulationInfo *sinfo) {
  const double dt2 = sinfo->TimeStep * 0.5;
  const int pn = vars->GetParticleNumber();
  double (*p)[D] = vars->p;
  const double exp1 = exp(-dt2 * vars->Zeta);
  for (int i = 0; i < pn; i++) {
    for (int d = 0; d < D; d++) {
      p[i][d] *= exp1;
    }
  }
}
//----------------------------------------------------------------------
void
ForceCalculator::Langevin(Variables *vars, SimulationInfo *sinfo, std::mt19937 &mt) {
  const double dt = sinfo->TimeStep;
  const int pn = vars->GetParticleNumber();
  double (*p)[D] = vars->p;
  const double hb_gamma = sinfo->HeatbathGamma;
  const double T = sinfo->AimedTemperature;
  const double hb_D = sqrt(2.0 * hb_gamma * T / dt);
  std::normal_distribution<double> nd(0.0,hb_D);
  for (int i = 0; i < pn; i++) {
    p[i][X] += (-hb_gamma * p[i][X] + nd(mt)) * dt;
    p[i][Y] += (-hb_gamma * p[i][Y] + nd(mt)) * dt;
    p[i][Z] += (-hb_gamma * p[i][Z] + nd(mt)) * dt;
  }
}
//----------------------------------------------------------------------
