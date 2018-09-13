//----------------------------------------------------------------------
#ifndef simd_h
#define simd_h
#include <emmintrin.h>
//----------------------------------------------------------------------
class v2df {
public:
  __m128d value;
  v2df(const double &v1, const double &v2) {value = _mm_set_pd(v2, v1);}
  v2df(const __m128d& v) {value = v;}
  v2df(const v2df& v0) {value = v0.value;}
  void set_val(const double &v1, const double &v2) {value = _mm_set_pd(v2, v1);}
  v2df operator+(const v2df& v0) {return v2df(__builtin_fj_add_v2r8(value, v0.value));}
  v2df operator-(const v2df& v0) {return v2df(__builtin_fj_sub_v2r8(value, v0.value));}
  v2df operator*(const v2df& v0) {return v2df(__builtin_fj_mul_v2r8(value, v0.value));}
  v2df operator/(const v2df& v0) {return v2df(__builtin_fj_mul_v2r8(value, __builtin_fj_rcpa_v2r8(v0.value)));}
  v2df inv() {return v2df(__builtin_fj_rcpa_v2r8(value));}
  v2df operator+=(const v2df& v0) {value = __builtin_fj_add_v2r8(value, v0.value); return *this;}
  v2df operator-=(const v2df& v0) {value = __builtin_fj_sub_v2r8(value, v0.value); return *this;}
  v2df max(const v2df& v0)       {return v2df(__builtin_fj_max_v2r8(value, v0.value));}
  v2df min(const v2df& v0)       {return v2df(__builtin_fj_min_v2r8(value, v0.value));}
  v2df madd(const v2df& v1, v2df& v2) {return v2df(__builtin_fj_madd_v2r8(value, v1.value, v2.value));}
  v2df msub(const v2df& v1, v2df& v2) {return v2df(__builtin_fj_msub_v2r8(value, v1.value, v2.value));}
  v2df floor() {return v2df(_fjsp_xtod_v2r8(_fjsp_dtox_v2r8(value)));}
  double operator[] (const int i) {double *p = (double*)(&value); return p[i];}
  double* pt() {return (double*)(&value);}
};
//----------------------------------------------------------------------
#endif
