//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: qscmvnvR.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 14-Apr-2020 21:40:05
//

// Include Files
#include "qscmvnvR.h"
#include "diag.h"
#include "eig.h"
#include "eml_primes_core.h"
#include "mtimes.h"
#include "qscmvnvR_data.h"
#include "qscmvnvR_emxutil.h"
#include "qscmvnvR_initialize.h"
#include "qscmvnvR_rtwutil.h"
#include "rand.h"
#include "relop.h"
#include "rt_nonfinite.h"
#include "sort.h"
#include "sqrt.h"
#include "sum.h"
#include <cmath>
#include <cstring>
#include <math.h>

// Function Declarations
static void chlsrt(const emxArray_real_T *r, emxArray_real_T *a, const
                   emxArray_real_T *cn, emxArray_real_T *b, emxArray_creal_T *ch,
                   emxArray_real_T *clg, double *np);
static double eml_erfcore(double x);
static void mvndnv(double n, const emxArray_real_T *a, const emxArray_creal_T
                   *ch, const emxArray_real_T *b, const emxArray_real_T *clg,
                   double ci, double dci, const emxArray_real_T *x, double nv,
                   emxArray_real_T *p);
static double rt_powd_snf(double u0, double u1);

// Function Definitions

//
// Computes permuted lower Cholesky factor ch for covariance r which
//    may be singular, combined with contraints a < cn*x < b, to
//    form revised lower triangular constraint set ap < ch*x < bp;
//    clg contains information about structure of ch: clg(1) rows for
//    ch with 1 nonzero, ..., clg(np) rows with np nonzeros.
// Arguments    : const emxArray_real_T *r
//                emxArray_real_T *a
//                const emxArray_real_T *cn
//                emxArray_real_T *b
//                emxArray_creal_T *ch
//                emxArray_real_T *clg
//                double *np
// Return Type  : void
//
static void chlsrt(const emxArray_real_T *r, emxArray_real_T *a, const
                   emxArray_real_T *cn, emxArray_real_T *b, emxArray_creal_T *ch,
                   emxArray_real_T *clg, double *np)
{
  emxArray_real_T *y;
  int i;
  int loop_ub;
  emxArray_real_T *b_ch;
  int m;
  emxArray_real_T *c;
  emxArray_real_T *d;
  emxArray_real_T *varargin_1;
  int nx;
  int k;
  emxArray_real_T *b_c;
  int b_i;
  emxArray_creal_T *v;
  int i1;
  emxArray_creal_T *vt;
  emxArray_creal_T *V;
  emxArray_creal_T *D;
  emxArray_creal_T *b_d;
  emxArray_int32_T *iidx;
  emxArray_boolean_T *x;
  int nz;
  int b_y;
  int b_loop_ub;
  emxArray_creal_T *c_ch;
  emxArray_creal_T *b_V;
  double dsb;
  double dsa;
  int i2;
  emxArray_creal_T *c_y;
  emxArray_creal_T *d_ch;
  double epi;
  double vm;
  unsigned int jm;
  unsigned int lm;
  int l;
  int c_loop_ub;
  unsigned int b_l;
  double dna;
  int i3;
  int i4;
  int loop_ub_tmp;
  unsigned int unnamed_idx_1;
  int i5;
  creal_T ss;
  boolean_T exitg1;
  double mn;
  double brm;
  double a_re;
  double dnb;
  double b_ss;
  int i6;
  int i7;
  int i8;
  int i9;
  emxInit_real_T(&y, 1);

  //
  //  end mvndnv
  //
  //  singularity tolerance;
  //
  i = y->size[0];
  y->size[0] = r->size[1];
  emxEnsureCapacity_real_T(y, i);
  loop_ub = r->size[1];
  for (i = 0; i < loop_ub; i++) {
    y->data[i] = 0.0;
  }

  i = clg->size[0] * clg->size[1];
  clg->size[0] = 1;
  clg->size[1] = r->size[1];
  emxEnsureCapacity_real_T(clg, i);
  loop_ub = r->size[1];
  for (i = 0; i < loop_ub; i++) {
    clg->data[i] = 0.0;
  }

  emxInit_real_T(&b_ch, 2);
  m = cn->size[0] - 1;
  i = b_ch->size[0] * b_ch->size[1];
  b_ch->size[0] = cn->size[0];
  b_ch->size[1] = cn->size[1];
  emxEnsureCapacity_real_T(b_ch, i);
  loop_ub = cn->size[0] * cn->size[1];
  for (i = 0; i < loop_ub; i++) {
    b_ch->data[i] = cn->data[i];
  }

  emxInit_real_T(&c, 2);
  i = c->size[0] * c->size[1];
  c->size[0] = r->size[0];
  c->size[1] = r->size[1];
  emxEnsureCapacity_real_T(c, i);
  loop_ub = r->size[0] * r->size[1];
  for (i = 0; i < loop_ub; i++) {
    c->data[i] = r->data[i];
  }

  emxInit_real_T(&d, 1);
  emxInit_real_T(&varargin_1, 1);
  diag(r, varargin_1);
  i = d->size[0];
  d->size[0] = varargin_1->size[0];
  emxEnsureCapacity_real_T(d, i);
  nx = varargin_1->size[0];
  for (k = 0; k < nx; k++) {
    if (varargin_1->data[k] > 0.0) {
      d->data[k] = varargin_1->data[k];
    } else {
      d->data[k] = 0.0;
    }
  }

  nx = d->size[0];
  for (k = 0; k < nx; k++) {
    d->data[k] = std::sqrt(d->data[k]);
  }

  i = cn->size[1];
  emxInit_real_T(&b_c, 2);
  for (b_i = 0; b_i < i; b_i++) {
    if (d->data[b_i] > 0.0) {
      loop_ub = c->size[0] - 1;
      i1 = varargin_1->size[0];
      varargin_1->size[0] = c->size[0];
      emxEnsureCapacity_real_T(varargin_1, i1);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        varargin_1->data[i1] = c->data[i1 + c->size[0] * b_i] / d->data[b_i];
      }

      loop_ub = varargin_1->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        c->data[i1 + c->size[0] * b_i] = varargin_1->data[i1];
      }

      loop_ub = c->size[1] - 1;
      i1 = b_c->size[0] * b_c->size[1];
      b_c->size[0] = 1;
      b_c->size[1] = c->size[1];
      emxEnsureCapacity_real_T(b_c, i1);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        b_c->data[i1] = c->data[b_i + c->size[0] * i1] / d->data[b_i];
      }

      loop_ub = b_c->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        c->data[b_i + c->size[0] * i1] = b_c->data[i1];
      }

      loop_ub = b_ch->size[0] - 1;
      i1 = varargin_1->size[0];
      varargin_1->size[0] = b_ch->size[0];
      emxEnsureCapacity_real_T(varargin_1, i1);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        varargin_1->data[i1] = b_ch->data[i1 + b_ch->size[0] * b_i] * d->
          data[b_i];
      }

      loop_ub = varargin_1->size[0];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_ch->data[i1 + b_ch->size[0] * b_i] = varargin_1->data[i1];
      }
    }
  }

  emxFree_real_T(&varargin_1);
  emxFree_real_T(&d);
  emxInit_creal_T(&v, 2);
  emxInit_creal_T(&vt, 1);
  emxInit_creal_T(&V, 2);
  emxInit_creal_T(&D, 2);

  //
  //   determine r factors and form revised constraint matrix ch
  //
  eig(c, V, D);
  b_diag(D, vt);
  i = v->size[0] * v->size[1];
  v->size[0] = 1;
  v->size[1] = vt->size[0];
  emxEnsureCapacity_creal_T(v, i);
  loop_ub = vt->size[0];
  emxFree_real_T(&c);
  for (i = 0; i < loop_ub; i++) {
    v->data[i].re = vt->data[i].re;
    v->data[i].im = -vt->data[i].im;
  }

  emxInit_creal_T(&b_d, 2);
  emxInit_int32_T(&iidx, 2);
  sort(v, iidx);
  i = b_d->size[0] * b_d->size[1];
  b_d->size[0] = 1;
  b_d->size[1] = v->size[1];
  emxEnsureCapacity_creal_T(b_d, i);
  nx = v->size[1];
  for (k = 0; k < nx; k++) {
    b_d->data[k] = v->data[k];
    if (relop(v->data[k], 0.0)) {
      b_d->data[k].re = 0.0;
      b_d->data[k].im = 0.0;
    }
  }

  nx = b_d->size[1];
  for (k = 0; k < nx; k++) {
    b_sqrt(&b_d->data[k]);
  }

  emxInit_boolean_T(&x, 2);
  i = x->size[0] * x->size[1];
  x->size[0] = 1;
  x->size[1] = b_d->size[1];
  emxEnsureCapacity_boolean_T(x, i);
  loop_ub = b_d->size[0] * b_d->size[1];
  for (i = 0; i < loop_ub; i++) {
    x->data[i] = (b_d->data[i].re > 0.0);
  }

  nx = x->size[1];
  if (x->size[1] == 0) {
    nz = -1;
  } else {
    nz = x->data[0] - 1;
    for (k = 2; k <= nx; k++) {
      nz += x->data[k - 1];
    }
  }

  nx = x->size[1];
  if (x->size[1] == 0) {
    b_y = 0;
  } else {
    b_y = x->data[0];
    for (k = 2; k <= nx; k++) {
      b_y += x->data[k - 1];
    }
  }

  emxFree_boolean_T(&x);
  if (1 > b_y) {
    loop_ub = -1;
    b_loop_ub = -1;
  } else {
    loop_ub = nz;
    b_loop_ub = nz;
  }

  i = vt->size[0];
  vt->size[0] = cn->size[1];
  emxEnsureCapacity_creal_T(vt, i);
  k = cn->size[1];
  for (i = 0; i < k; i++) {
    vt->data[i].re = 1.0;
    vt->data[i].im = 0.0;
  }

  i = D->size[0] * D->size[1];
  D->size[0] = vt->size[0];
  D->size[1] = b_loop_ub + 1;
  emxEnsureCapacity_creal_T(D, i);
  for (i = 0; i <= b_loop_ub; i++) {
    k = vt->size[0];
    for (i1 = 0; i1 < k; i1++) {
      D->data[i1 + D->size[0] * i].re = vt->data[i1].re * b_d->data[i].re -
        vt->data[i1].im * b_d->data[i].im;
      D->data[i1 + D->size[0] * i].im = vt->data[i1].re * b_d->data[i].im +
        vt->data[i1].im * b_d->data[i].re;
    }
  }

  emxInit_creal_T(&c_ch, 2);
  i = c_ch->size[0] * c_ch->size[1];
  c_ch->size[0] = b_ch->size[0];
  c_ch->size[1] = b_ch->size[1];
  emxEnsureCapacity_creal_T(c_ch, i);
  b_loop_ub = b_ch->size[0] * b_ch->size[1];
  for (i = 0; i < b_loop_ub; i++) {
    c_ch->data[i].re = b_ch->data[i];
    c_ch->data[i].im = 0.0;
  }

  emxFree_real_T(&b_ch);
  emxInit_creal_T(&b_V, 2);
  b_loop_ub = V->size[0];
  i = b_V->size[0] * b_V->size[1];
  b_V->size[0] = V->size[0];
  b_V->size[1] = loop_ub + 1;
  emxEnsureCapacity_creal_T(b_V, i);
  for (i = 0; i <= loop_ub; i++) {
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      b_V->data[i1 + b_V->size[0] * i].re = V->data[i1 + V->size[0] *
        (iidx->data[i] - 1)].re * D->data[i1 + D->size[0] * i].re - V->data[i1 +
        V->size[0] * (iidx->data[i] - 1)].im * D->data[i1 + D->size[0] * i].im;
      b_V->data[i1 + b_V->size[0] * i].im = V->data[i1 + V->size[0] *
        (iidx->data[i] - 1)].re * D->data[i1 + D->size[0] * i].im + V->data[i1 +
        V->size[0] * (iidx->data[i] - 1)].im * D->data[i1 + D->size[0] * i].re;
    }
  }

  emxFree_int32_T(&iidx);
  emxFree_creal_T(&D);
  emxFree_creal_T(&V);
  i = ch->size[0] * ch->size[1];
  ch->size[0] = c_ch->size[0];
  ch->size[1] = b_V->size[1];
  emxEnsureCapacity_creal_T(ch, i);
  loop_ub = c_ch->size[0];
  for (i = 0; i < loop_ub; i++) {
    b_loop_ub = b_V->size[1];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      ch->data[i + ch->size[0] * i1].re = 0.0;
      ch->data[i + ch->size[0] * i1].im = 0.0;
      k = c_ch->size[1];
      for (i2 = 0; i2 < k; i2++) {
        ch->data[i + ch->size[0] * i1].re += c_ch->data[i + c_ch->size[0] * i2].
          re * b_V->data[i2 + b_V->size[0] * i1].re - c_ch->data[i + c_ch->size
          [0] * i2].im * b_V->data[i2 + b_V->size[0] * i1].im;
        ch->data[i + ch->size[0] * i1].im += c_ch->data[i + c_ch->size[0] * i2].
          re * b_V->data[i2 + b_V->size[0] * i1].im + c_ch->data[i + c_ch->size
          [0] * i2].im * b_V->data[i2 + b_V->size[0] * i1].re;
      }
    }
  }

  emxFree_creal_T(&b_V);

  //
  //  use right reflectors to reduce ch to lower triangular
  //
  dsb = static_cast<double>(b_y) - 1.0;
  dsa = cn->size[0];
  if (dsb < dsa) {
    dsa = dsb;
  }

  i = static_cast<int>(dsa);
  emxInit_creal_T(&c_y, 1);
  emxInit_creal_T(&d_ch, 2);
  for (b_i = 0; b_i < i; b_i++) {
    epi = 1.0E-10 * (static_cast<double>(b_i) + 1.0);
    vm = 1.0;
    jm = b_i + 1U;
    lm = jm;

    //
    //  permute rows so that smallest variance variables are first.
    //
    i1 = m - b_i;
    if (0 <= i1) {
      if (1 > b_i) {
        c_loop_ub = 0;
      } else {
        c_loop_ub = b_i;
      }

      if (static_cast<double>(b_i) + 1.0 > b_y) {
        i3 = 0;
        i4 = 0;
      } else {
        i3 = b_i;
        i4 = nz + 1;
      }

      loop_ub_tmp = i4 - i3;
      unnamed_idx_1 = static_cast<unsigned int>(loop_ub_tmp);
    }

    for (l = 0; l <= i1; l++) {
      b_l = (static_cast<unsigned int>(b_i) + l) + 1U;
      dna = 0.0;
      dsb = 0.0;
      for (i2 = 0; i2 < c_loop_ub; i2++) {
        k = static_cast<int>(b_l) - 1;
        dna += ch->data[k + ch->size[0] * i2].re * y->data[i2] - ch->data[k +
          ch->size[0] * i2].im * 0.0;
        dsb += ch->data[k + ch->size[0] * i2].re * 0.0 + ch->data[k + ch->size[0]
          * i2].im * y->data[i2];
      }

      i2 = v->size[0] * v->size[1];
      v->size[0] = 1;
      v->size[1] = i4 - i3;
      emxEnsureCapacity_creal_T(v, i2);
      for (i2 = 0; i2 < loop_ub_tmp; i2++) {
        v->data[i2] = ch->data[(static_cast<int>(b_l) + ch->size[0] * (i3 + i2))
          - 1];
      }

      i2 = b_d->size[0] * b_d->size[1];
      b_d->size[0] = 1;
      i5 = static_cast<int>(unnamed_idx_1);
      b_d->size[1] = static_cast<int>(unnamed_idx_1);
      emxEnsureCapacity_creal_T(b_d, i2);
      for (k = 0; k < i5; k++) {
        b_d->data[k].re = v->data[k].re * v->data[k].re - v->data[k].im *
          v->data[k].im;
        dsa = v->data[k].re * v->data[k].im;
        b_d->data[k].im = dsa + dsa;
      }

      ss = sum(b_d);
      b_sqrt(&ss);
      if (relop(ss, epi)) {
        ss.re = epi;
        ss.im = 0.0;
      }

      nx = static_cast<int>(b_l) - 1;
      dsa = a->data[nx] - dna;
      if (ss.im == 0.0) {
        if (0.0 - dsb == 0.0) {
          a_re = dsa / ss.re;
        } else if (dsa == 0.0) {
          a_re = 0.0;
        } else {
          a_re = dsa / ss.re;
        }
      } else if (ss.re == 0.0) {
        if (dsa == 0.0) {
          a_re = (0.0 - dsb) / ss.im;
        } else if (0.0 - dsb == 0.0) {
          a_re = 0.0;
        } else {
          a_re = (0.0 - dsb) / ss.im;
        }
      } else {
        brm = std::abs(ss.re);
        dnb = std::abs(ss.im);
        if (brm > dnb) {
          dnb = ss.im / ss.re;
          a_re = (dsa + dnb * (0.0 - dsb)) / (ss.re + dnb * ss.im);
        } else if (dnb == brm) {
          if (ss.re > 0.0) {
            b_ss = 0.5;
          } else {
            b_ss = -0.5;
          }

          if (ss.im > 0.0) {
            dnb = 0.5;
          } else {
            dnb = -0.5;
          }

          a_re = (dsa * b_ss + (0.0 - dsb) * dnb) / brm;
        } else {
          dnb = ss.re / ss.im;
          a_re = (dnb * dsa + (0.0 - dsb)) / (ss.im + dnb * ss.re);
        }
      }

      dsa = b->data[nx] - dna;
      if (ss.im == 0.0) {
        if (0.0 - dsb == 0.0) {
          brm = dsa / ss.re;
        } else if (dsa == 0.0) {
          brm = 0.0;
        } else {
          brm = dsa / ss.re;
        }
      } else if (ss.re == 0.0) {
        if (dsa == 0.0) {
          brm = (0.0 - dsb) / ss.im;
        } else if (0.0 - dsb == 0.0) {
          brm = 0.0;
        } else {
          brm = (0.0 - dsb) / ss.im;
        }
      } else {
        brm = std::abs(ss.re);
        dnb = std::abs(ss.im);
        if (brm > dnb) {
          dnb = ss.im / ss.re;
          brm = (dsa + dnb * (0.0 - dsb)) / (ss.re + dnb * ss.im);
        } else if (dnb == brm) {
          if (ss.re > 0.0) {
            b_ss = 0.5;
          } else {
            b_ss = -0.5;
          }

          if (ss.im > 0.0) {
            dnb = 0.5;
          } else {
            dnb = -0.5;
          }

          brm = (dsa * b_ss + (0.0 - dsb) * dnb) / brm;
        } else {
          dnb = ss.re / ss.im;
          brm = (dnb * dsa + (0.0 - dsb)) / (ss.im + dnb * ss.re);
        }
      }

      dna = 0.0;
      dsa = 0.0;
      dnb = 0.0;
      dsb = 1.0;
      if (a_re > -9.0) {
        dna = std::exp(-0.5 * a_re * a_re) / 2.5066282746310002;
        mn = eml_erfcore(-a_re / 1.4142135623730951);
        dsa = 0.5 * mn;
      }

      if (brm < 9.0) {
        dnb = std::exp(-0.5 * brm * brm) / 2.5066282746310002;
        mn = eml_erfcore(-brm / 1.4142135623730951);
        dsb = 0.5 * mn;
      }

      dsa = dsb - dsa;
      if (dsa > epi) {
        if (a_re <= -9.0) {
          mn = -dnb;
          dsb = -brm * dnb;
        } else if (brm >= 9.0) {
          mn = dna;
          dsb = a_re * dna;
        } else {
          mn = dna - dnb;
          dsb = a_re * dna - brm * dnb;
        }

        mn /= dsa;
        dsb = (dsb / dsa + 1.0) - mn * mn;
      } else {
        mn = (a_re + brm) / 2.0;
        dsb = 0.0;
        if (a_re <= -9.0) {
          mn = brm;
        } else {
          if (brm >= 9.0) {
            mn = a_re;
          }
        }
      }

      if (dsb <= vm) {
        lm = b_l;
        vm = dsb;
        y->data[b_i] = mn;
      }
    }

    if (1 > b_y) {
      loop_ub = -1;
    } else {
      loop_ub = nz;
    }

    i1 = v->size[0] * v->size[1];
    v->size[0] = 1;
    v->size[1] = loop_ub + 1;
    emxEnsureCapacity_creal_T(v, i1);
    for (i1 = 0; i1 <= loop_ub; i1++) {
      v->data[i1] = ch->data[(static_cast<int>(lm) + ch->size[0] * i1) - 1];
    }

    if (lm > jm) {
      if (1 > b_y) {
        loop_ub = -1;
      } else {
        loop_ub = nz;
      }

      nx = static_cast<int>(lm) - 1;
      i1 = d_ch->size[0] * d_ch->size[1];
      d_ch->size[0] = 2;
      d_ch->size[1] = loop_ub + 1;
      emxEnsureCapacity_creal_T(d_ch, i1);
      for (i1 = 0; i1 <= loop_ub; i1++) {
        d_ch->data[2 * i1] = ch->data[nx + ch->size[0] * i1];
        d_ch->data[2 * i1 + 1] = ch->data[b_i + ch->size[0] * i1];
      }

      loop_ub = d_ch->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        ch->data[b_i + ch->size[0] * i1] = d_ch->data[2 * i1];
        ch->data[nx + ch->size[0] * i1] = d_ch->data[2 * i1 + 1];
      }

      dsb = a->data[b_i];
      a->data[b_i] = a->data[nx];
      a->data[nx] = dsb;
      dsb = b->data[b_i];
      b->data[b_i] = b->data[nx];
      b->data[nx] = dsb;
    }

    if ((static_cast<double>(b_i) + 1.0) + 1.0 > b_y) {
      i1 = 0;
      i2 = 0;
    } else {
      i1 = b_i + 1;
      i2 = nz + 1;
    }

    loop_ub = i2 - i1;
    for (i2 = 0; i2 < loop_ub; i2++) {
      i5 = i1 + i2;
      ch->data[b_i + ch->size[0] * i5].re = 0.0;
      ch->data[b_i + ch->size[0] * i5].im = 0.0;
    }

    if ((static_cast<double>(b_i) + 1.0) + 1.0 > b_y) {
      i1 = -1;
      i2 = -1;
    } else {
      i1 = b_i;
      i2 = nz;
    }

    unnamed_idx_1 = static_cast<unsigned int>((i2 - i1));
    i2 = b_d->size[0] * b_d->size[1];
    b_d->size[0] = 1;
    b_d->size[1] = static_cast<int>(unnamed_idx_1);
    emxEnsureCapacity_creal_T(b_d, i2);
    nx = static_cast<int>(unnamed_idx_1);
    for (k = 0; k < nx; k++) {
      i2 = (i1 + k) + 1;
      b_d->data[k].re = v->data[i2].re * v->data[i2].re - v->data[i2].im *
        v->data[i2].im;
      dsa = v->data[i2].re * v->data[i2].im;
      b_d->data[k].im = dsa + dsa;
    }

    ss = sum(b_d);
    if (ss.re > epi) {
      dsb = v->data[b_i].re * v->data[b_i].im;
      ss.re += v->data[b_i].re * v->data[b_i].re - v->data[b_i].im * v->data[b_i]
        .im;
      ss.im += dsb + dsb;
      b_sqrt(&ss);
      if (v->data[b_i].re < 0.0) {
        ss.re = -ss.re;
        ss.im = -ss.im;
      }

      ch->data[b_i + ch->size[0] * b_i].re = -ss.re;
      ch->data[b_i + ch->size[0] * b_i].im = -ss.im;
      v->data[b_i].re += ss.re;
      v->data[b_i].im += ss.im;
      if (static_cast<double>(b_i) + 1.0 > b_y) {
        i1 = 0;
        i2 = -1;
      } else {
        i1 = b_i;
        i2 = nz;
      }

      mn = ss.re * v->data[b_i].re - ss.im * v->data[b_i].im;
      dna = ss.re * v->data[b_i].im + ss.im * v->data[b_i].re;
      loop_ub = i2 - i1;
      i2 = vt->size[0];
      vt->size[0] = loop_ub + 1;
      emxEnsureCapacity_creal_T(vt, i2);
      for (i2 = 0; i2 <= loop_ub; i2++) {
        nx = i1 + i2;
        if (dna == 0.0) {
          if (-v->data[nx].im == 0.0) {
            vt->data[i2].re = v->data[nx].re / mn;
            vt->data[i2].im = 0.0;
          } else if (v->data[nx].re == 0.0) {
            vt->data[i2].re = 0.0;
            vt->data[i2].im = -v->data[nx].im / mn;
          } else {
            vt->data[i2].re = v->data[nx].re / mn;
            vt->data[i2].im = -v->data[nx].im / mn;
          }
        } else if (mn == 0.0) {
          if (v->data[nx].re == 0.0) {
            vt->data[i2].re = -v->data[nx].im / dna;
            vt->data[i2].im = 0.0;
          } else if (-v->data[nx].im == 0.0) {
            vt->data[i2].re = 0.0;
            vt->data[i2].im = -(v->data[nx].re / dna);
          } else {
            vt->data[i2].re = -v->data[nx].im / dna;
            vt->data[i2].im = -(v->data[nx].re / dna);
          }
        } else {
          brm = std::abs(mn);
          dnb = std::abs(dna);
          if (brm > dnb) {
            dnb = dna / mn;
            dsb = mn + dnb * dna;
            vt->data[i2].re = (v->data[nx].re + dnb * -v->data[nx].im) / dsb;
            vt->data[i2].im = (-v->data[nx].im - dnb * v->data[nx].re) / dsb;
          } else if (dnb == brm) {
            if (mn > 0.0) {
              dsa = 0.5;
            } else {
              dsa = -0.5;
            }

            if (dna > 0.0) {
              dsb = 0.5;
            } else {
              dsb = -0.5;
            }

            vt->data[i2].re = (v->data[nx].re * dsa + -v->data[nx].im * dsb) /
              brm;
            vt->data[i2].im = (-v->data[nx].im * dsa - v->data[nx].re * dsb) /
              brm;
          } else {
            dnb = mn / dna;
            dsb = dna + dnb * mn;
            vt->data[i2].re = (dnb * v->data[nx].re + -v->data[nx].im) / dsb;
            vt->data[i2].im = (dnb * -v->data[nx].im - v->data[nx].re) / dsb;
          }
        }
      }

      jm = b_i + 2U;
      if (jm > static_cast<unsigned int>((m + 1))) {
        i1 = -1;
      } else {
        i1 = b_i;
      }

      if (static_cast<double>(b_i) + 1.0 > b_y) {
        i2 = 0;
      } else {
        i2 = b_i;
      }

      if (jm > static_cast<unsigned int>((m + 1))) {
        i5 = 0;
        l = -1;
      } else {
        i5 = b_i + 1;
        l = m;
      }

      if (static_cast<double>(b_i) + 1.0 > b_y) {
        i6 = 0;
        i7 = -1;
        i8 = 0;
        i9 = -1;
      } else {
        i6 = b_i;
        i7 = nz;
        i8 = b_i;
        i9 = nz;
      }

      loop_ub = i7 - i6;
      i7 = loop_ub + 1;
      if ((i7 == 1) || (vt->size[0] == 1)) {
        b_loop_ub = l - i5;
        l = c_y->size[0];
        c_y->size[0] = b_loop_ub + 1;
        emxEnsureCapacity_creal_T(c_y, l);
        for (l = 0; l <= b_loop_ub; l++) {
          c_y->data[l].re = 0.0;
          c_y->data[l].im = 0.0;
          for (i7 = 0; i7 <= loop_ub; i7++) {
            k = i5 + l;
            nx = i6 + i7;
            c_y->data[l].re += ch->data[k + ch->size[0] * nx].re * vt->data[i7].
              re - ch->data[k + ch->size[0] * nx].im * vt->data[i7].im;
            c_y->data[l].im += ch->data[k + ch->size[0] * nx].re * vt->data[i7].
              im + ch->data[k + ch->size[0] * nx].im * vt->data[i7].re;
          }
        }
      } else {
        nx = c_ch->size[0] * c_ch->size[1];
        b_loop_ub = l - i5;
        c_ch->size[0] = b_loop_ub + 1;
        c_ch->size[1] = i7;
        emxEnsureCapacity_creal_T(c_ch, nx);
        for (l = 0; l <= loop_ub; l++) {
          for (i7 = 0; i7 <= b_loop_ub; i7++) {
            c_ch->data[i7 + c_ch->size[0] * l] = ch->data[(i5 + i7) + ch->size[0]
              * (i6 + l)];
          }
        }

        mtimes(c_ch, vt, c_y);
      }

      if (jm > static_cast<unsigned int>((m + 1))) {
        i5 = -1;
      } else {
        i5 = b_i;
      }

      if (static_cast<double>(b_i) + 1.0 > b_y) {
        l = 0;
      } else {
        l = b_i;
      }

      i6 = c_ch->size[0] * c_ch->size[1];
      c_ch->size[0] = c_y->size[0];
      loop_ub = i9 - i8;
      c_ch->size[1] = loop_ub + 1;
      emxEnsureCapacity_creal_T(c_ch, i6);
      b_loop_ub = c_y->size[0];
      for (i6 = 0; i6 < b_loop_ub; i6++) {
        for (i7 = 0; i7 <= loop_ub; i7++) {
          k = i8 + i7;
          i9 = (i1 + i6) + 1;
          nx = i2 + i7;
          c_ch->data[i6 + c_ch->size[0] * i7].re = ch->data[i9 + ch->size[0] *
            nx].re - (c_y->data[i6].re * v->data[k].re - c_y->data[i6].im *
                      v->data[k].im);
          c_ch->data[i6 + c_ch->size[0] * i7].im = ch->data[i9 + ch->size[0] *
            nx].im - (c_y->data[i6].re * v->data[k].im + c_y->data[i6].im *
                      v->data[k].re);
        }
      }

      loop_ub = c_ch->size[1];
      for (i1 = 0; i1 < loop_ub; i1++) {
        b_loop_ub = c_ch->size[0];
        for (i2 = 0; i2 < b_loop_ub; i2++) {
          ch->data[((i5 + i2) + ch->size[0] * (l + i1)) + 1] = c_ch->data[i2 +
            c_ch->size[0] * i1];
        }
      }
    }
  }

  emxFree_creal_T(&d_ch);
  emxFree_creal_T(&c_y);
  emxFree_creal_T(&vt);

  //
  //  scale and sort constraints
  //
  i = y->size[0];
  y->size[0] = cn->size[0];
  emxEnsureCapacity_real_T(y, i);
  loop_ub = cn->size[0];
  for (i = 0; i < loop_ub; i++) {
    y->data[i] = 0.0;
  }

  //  added
  i = cn->size[0];
  for (b_i = 0; b_i < i; b_i++) {
    if (1 > b_y) {
      loop_ub = 0;
    } else {
      loop_ub = nz + 1;
    }

    i1 = v->size[0] * v->size[1];
    v->size[0] = 1;
    v->size[1] = loop_ub;
    emxEnsureCapacity_creal_T(v, i1);
    for (i1 = 0; i1 < loop_ub; i1++) {
      v->data[i1] = ch->data[b_i + ch->size[0] * i1];
    }

    if (b_i + 1 < nz + 1) {
      y->data[b_i] = static_cast<double>(b_i) + 1.0;
    } else {
      y->data[b_i] = nz + 1;
    }

    jm = 1U;
    i1 = static_cast<int>(y->data[b_i]);
    for (nx = 0; nx < i1; nx++) {
      if (rt_hypotd_snf(v->data[nx].re, v->data[nx].im) > 1.0E-10 * (
           static_cast<double>(nx) + 1.0)) {
        jm = nx + 1U;
      }
    }

    if (static_cast<double>(jm) < nz + 1) {
      if (static_cast<int>(jm) + 1 > b_y) {
        i1 = 0;
        i2 = -1;
      } else {
        i1 = static_cast<int>(jm);
        i2 = nz;
      }

      loop_ub = (i2 - i1) + 1;
      for (i2 = 0; i2 < loop_ub; i2++) {
        i3 = i1 + i2;
        v->data[i3].re = 0.0;
        v->data[i3].im = 0.0;
      }
    }

    k = static_cast<int>(jm) - 1;
    clg->data[k]++;
    mn = a->data[b_i];
    dna = b->data[b_i];
    nx = b_i;
    i1 = static_cast<int>((((-1.0 - ((static_cast<double>(b_i) + 1.0) - 1.0)) +
      1.0) / -1.0));
    l = 0;
    exitg1 = false;
    while ((!exitg1) && (l <= i1 - 1)) {
      i2 = b_i - l;
      i3 = i2 - 1;
      if (jm >= y->data[i3]) {
        exitg1 = true;
      } else {
        if (1 > b_y) {
          loop_ub = -1;
        } else {
          loop_ub = nz;
        }

        i4 = b_d->size[0] * b_d->size[1];
        b_d->size[0] = 1;
        b_d->size[1] = loop_ub + 1;
        emxEnsureCapacity_creal_T(b_d, i4);
        for (i4 = 0; i4 <= loop_ub; i4++) {
          b_d->data[i4] = ch->data[i3 + ch->size[0] * i4];
        }

        loop_ub = b_d->size[1];
        for (i4 = 0; i4 < loop_ub; i4++) {
          ch->data[i2 + ch->size[0] * i4] = b_d->data[i4];
        }

        nx = i3;
        a->data[i2] = a->data[i3];
        b->data[i2] = b->data[i3];
        y->data[i2] = y->data[i3];
        l++;
      }
    }

    y->data[nx] = jm;
    loop_ub = v->size[1];
    for (i1 = 0; i1 < loop_ub; i1++) {
      if (v->data[k].im == 0.0) {
        if (v->data[i1].im == 0.0) {
          ch->data[nx + ch->size[0] * i1].re = v->data[i1].re / v->data[k].re;
          ch->data[nx + ch->size[0] * i1].im = 0.0;
        } else if (v->data[i1].re == 0.0) {
          ch->data[nx + ch->size[0] * i1].re = 0.0;
          ch->data[nx + ch->size[0] * i1].im = v->data[i1].im / v->data[k].re;
        } else {
          ch->data[nx + ch->size[0] * i1].re = v->data[i1].re / v->data[k].re;
          ch->data[nx + ch->size[0] * i1].im = v->data[i1].im / v->data[k].re;
        }
      } else if (v->data[k].re == 0.0) {
        if (v->data[i1].re == 0.0) {
          ch->data[nx + ch->size[0] * i1].re = v->data[i1].im / v->data[k].im;
          ch->data[nx + ch->size[0] * i1].im = 0.0;
        } else if (v->data[i1].im == 0.0) {
          ch->data[nx + ch->size[0] * i1].re = 0.0;
          ch->data[nx + ch->size[0] * i1].im = -(v->data[i1].re / v->data[k].im);
        } else {
          ch->data[nx + ch->size[0] * i1].re = v->data[i1].im / v->data[k].im;
          ch->data[nx + ch->size[0] * i1].im = -(v->data[i1].re / v->data[k].im);
        }
      } else {
        brm = std::abs(v->data[k].re);
        dnb = std::abs(v->data[k].im);
        if (brm > dnb) {
          dnb = v->data[k].im / v->data[k].re;
          dsb = v->data[k].re + dnb * v->data[k].im;
          ch->data[nx + ch->size[0] * i1].re = (v->data[i1].re + dnb * v->
            data[i1].im) / dsb;
          ch->data[nx + ch->size[0] * i1].im = (v->data[i1].im - dnb * v->
            data[i1].re) / dsb;
        } else if (dnb == brm) {
          if (v->data[k].re > 0.0) {
            dsa = 0.5;
          } else {
            dsa = -0.5;
          }

          if (v->data[k].im > 0.0) {
            dsb = 0.5;
          } else {
            dsb = -0.5;
          }

          ch->data[nx + ch->size[0] * i1].re = (v->data[i1].re * dsa + v->
            data[i1].im * dsb) / brm;
          ch->data[nx + ch->size[0] * i1].im = (v->data[i1].im * dsa - v->
            data[i1].re * dsb) / brm;
        } else {
          dnb = v->data[k].re / v->data[k].im;
          dsb = v->data[k].im + dnb * v->data[k].re;
          ch->data[nx + ch->size[0] * i1].re = (dnb * v->data[i1].re + v->
            data[i1].im) / dsb;
          ch->data[nx + ch->size[0] * i1].im = (dnb * v->data[i1].im - v->
            data[i1].re) / dsb;
        }
      }
    }

    if (v->data[k].im == 0.0) {
      dsb = mn / v->data[k].re;
    } else if (v->data[k].re == 0.0) {
      if (mn == 0.0) {
        dsb = 0.0 / v->data[k].im;
      } else {
        dsb = 0.0;
      }
    } else {
      brm = std::abs(v->data[k].re);
      dnb = std::abs(v->data[k].im);
      if (brm > dnb) {
        dnb = v->data[k].im / v->data[k].re;
        dsb = (mn + dnb * 0.0) / (v->data[k].re + dnb * v->data[k].im);
      } else if (dnb == brm) {
        if (v->data[k].re > 0.0) {
          dnb = 0.5;
        } else {
          dnb = -0.5;
        }

        if (v->data[k].im > 0.0) {
          b_ss = 0.5;
        } else {
          b_ss = -0.5;
        }

        dsb = (mn * dnb + 0.0 * b_ss) / brm;
      } else {
        dnb = v->data[k].re / v->data[k].im;
        dsb = dnb * mn / (v->data[k].im + dnb * v->data[k].re);
      }
    }

    a->data[nx] = dsb;
    if (v->data[k].im == 0.0) {
      dsb = dna / v->data[k].re;
    } else if (v->data[k].re == 0.0) {
      if (dna == 0.0) {
        dsb = 0.0 / v->data[k].im;
      } else {
        dsb = 0.0;
      }
    } else {
      brm = std::abs(v->data[k].re);
      dnb = std::abs(v->data[k].im);
      if (brm > dnb) {
        dnb = v->data[k].im / v->data[k].re;
        dsb = (dna + dnb * 0.0) / (v->data[k].re + dnb * v->data[k].im);
      } else if (dnb == brm) {
        if (v->data[k].re > 0.0) {
          dnb = 0.5;
        } else {
          dnb = -0.5;
        }

        if (v->data[k].im > 0.0) {
          b_ss = 0.5;
        } else {
          b_ss = -0.5;
        }

        dsb = (dna * dnb + 0.0 * b_ss) / brm;
      } else {
        dnb = v->data[k].re / v->data[k].im;
        dsb = dnb * dna / (v->data[k].im + dnb * v->data[k].re);
      }
    }

    b->data[nx] = dsb;
    if (v->data[k].re < 0.0) {
      dsb = a->data[nx];
      a->data[nx] = b->data[nx];
      b->data[nx] = dsb;
    }
  }

  emxFree_creal_T(&b_d);
  emxFree_creal_T(&v);
  emxFree_real_T(&y);
  jm = 0U;
  for (b_i = 0; b_i <= nz; b_i++) {
    if (clg->data[b_i] > 0.0) {
      jm = b_i + 1U;
    }
  }

  //
  //  combine constraints for first variable
  //
  if (clg->data[0] > 1.0) {
    if (1.0 > clg->data[0]) {
      i = 0;
    } else {
      i = static_cast<int>(clg->data[0]);
    }

    if (i <= 2) {
      if ((i != 1) && ((a->data[0] < a->data[1]) || (rtIsNaN(a->data[0]) &&
            (!rtIsNaN(a->data[1]))))) {
        a->data[0] = a->data[1];
      }
    } else {
      if (!rtIsNaN(a->data[0])) {
        nx = 1;
      } else {
        nx = 0;
        k = 2;
        exitg1 = false;
        while ((!exitg1) && (k <= i)) {
          if (!rtIsNaN(a->data[k - 1])) {
            nx = k;
            exitg1 = true;
          } else {
            k++;
          }
        }
      }

      if (nx != 0) {
        dsb = a->data[nx - 1];
        i1 = nx + 1;
        for (k = i1; k <= i; k++) {
          dsa = a->data[k - 1];
          if (dsb < dsa) {
            dsb = dsa;
          }
        }

        a->data[0] = dsb;
      }
    }

    if (1.0 > clg->data[0]) {
      i = 0;
    } else {
      i = static_cast<int>(clg->data[0]);
    }

    if (i <= 2) {
      if (i == 1) {
        dsb = b->data[0];
      } else if ((b->data[0] > b->data[1]) || (rtIsNaN(b->data[0]) && (!rtIsNaN
                   (b->data[1])))) {
        dsb = b->data[1];
      } else {
        dsb = b->data[0];
      }
    } else {
      if (!rtIsNaN(b->data[0])) {
        nx = 1;
      } else {
        nx = 0;
        k = 2;
        exitg1 = false;
        while ((!exitg1) && (k <= i)) {
          if (!rtIsNaN(b->data[k - 1])) {
            nx = k;
            exitg1 = true;
          } else {
            k++;
          }
        }
      }

      if (nx == 0) {
        dsb = b->data[0];
      } else {
        dsb = b->data[nx - 1];
        i1 = nx + 1;
        for (k = i1; k <= i; k++) {
          dsa = b->data[k - 1];
          if (dsb > dsa) {
            dsb = dsa;
          }
        }
      }
    }

    if ((a->data[0] > dsb) || rtIsNaN(dsb)) {
      b->data[0] = a->data[0];
    } else {
      b->data[0] = dsb;
    }

    if (clg->data[0] + 1.0 > cn->size[0]) {
      i = 0;
      i1 = 0;
    } else {
      i = static_cast<int>((clg->data[0] + 1.0)) - 1;
      i1 = cn->size[0];
    }

    i2 = !(2.0 > (static_cast<double>(cn->size[0]) - clg->data[0]) + 1.0);
    i3 = b_c->size[0] * b_c->size[1];
    b_c->size[0] = 1;
    loop_ub = i1 - i;
    b_c->size[1] = loop_ub;
    emxEnsureCapacity_real_T(b_c, i3);
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_c->data[i1] = a->data[i + i1];
    }

    loop_ub = b_c->size[1];
    for (i = 0; i < loop_ub; i++) {
      a->data[i2 + i] = b_c->data[i];
    }

    if (clg->data[0] + 1.0 > cn->size[0]) {
      i = 0;
      i1 = 0;
    } else {
      i = static_cast<int>((clg->data[0] + 1.0)) - 1;
      i1 = cn->size[0];
    }

    i2 = !(2.0 > (static_cast<double>(cn->size[0]) - clg->data[0]) + 1.0);
    i3 = b_c->size[0] * b_c->size[1];
    b_c->size[0] = 1;
    loop_ub = i1 - i;
    b_c->size[1] = loop_ub;
    emxEnsureCapacity_real_T(b_c, i3);
    for (i1 = 0; i1 < loop_ub; i1++) {
      b_c->data[i1] = b->data[i + i1];
    }

    loop_ub = b_c->size[1];
    for (i = 0; i < loop_ub; i++) {
      b->data[i2 + i] = b_c->data[i];
    }

    if (clg->data[0] + 1.0 > cn->size[0]) {
      i = 0;
      i1 = 0;
    } else {
      i = static_cast<int>((clg->data[0] + 1.0)) - 1;
      i1 = cn->size[0];
    }

    i2 = !(2.0 > (static_cast<double>(cn->size[0]) - clg->data[0]) + 1.0);
    loop_ub = ch->size[1] - 1;
    b_loop_ub = i1 - i;
    i1 = c_ch->size[0] * c_ch->size[1];
    c_ch->size[0] = b_loop_ub;
    c_ch->size[1] = ch->size[1];
    emxEnsureCapacity_creal_T(c_ch, i1);
    for (i1 = 0; i1 <= loop_ub; i1++) {
      for (i3 = 0; i3 < b_loop_ub; i3++) {
        c_ch->data[i3 + c_ch->size[0] * i1] = ch->data[(i + i3) + ch->size[0] *
          i1];
      }
    }

    loop_ub = c_ch->size[1];
    for (i = 0; i < loop_ub; i++) {
      b_loop_ub = c_ch->size[0];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        ch->data[(i2 + i1) + ch->size[0] * i] = c_ch->data[i1 + c_ch->size[0] *
          i];
      }
    }

    clg->data[0] = 1.0;
  }

  emxFree_creal_T(&c_ch);
  emxFree_real_T(&b_c);

  //
  //  end chlsrt
  //
  //   Standard statistical normal distribution pdf, cdf, and inverse
  // function p = phi(z),   p = normpdf(z);         % p = exp(-z.^2/2)/sqrt(2*pi); 
  // function p = Phi(z),   p = normcdf(z);         % p = erfc( -z/sqrt(2) )/2;  
  // function z = Phinv(w), z = norminv(w);         % z = -sqrt(2)*erfcinv( 2*w ); 
  *np = jm;
}

//
// Arguments    : double x
// Return Type  : double
//
static double eml_erfcore(double x)
{
  double y;
  double absx;
  double S;
  double s;
  double R;
  int eint;

  // ========================== COPYRIGHT NOTICE ============================
  //  The algorithms for calculating ERF(X) and ERFC(X) are derived
  //  from FDLIBM, which has the following notice:
  //
  //  Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.
  //
  //  Developed at SunSoft, a Sun Microsystems, Inc. business.
  //  Permission to use, copy, modify, and distribute this
  //  software is freely granted, provided that this notice
  //  is preserved.
  // =============================    END    ================================
  absx = std::abs(x);
  if (rtIsNaN(x)) {
    y = x;
  } else if (rtIsInf(x)) {
    if (x < 0.0) {
      y = 2.0;
    } else {
      y = 0.0;
    }
  } else if (absx < 0.84375) {
    if (absx < 1.3877787807814457E-17) {
      y = 1.0 - x;
    } else {
      s = x * x;
      y = (s * (s * (s * (s * -2.3763016656650163E-5 + -0.0057702702964894416) +
                     -0.02848174957559851) + -0.3250421072470015) +
           0.12837916709551256) / (s * (s * (s * (s * (s *
        -3.9602282787753681E-6 + 0.00013249473800432164) + 0.0050813062818757656)
        + 0.0650222499887673) + 0.39791722395915535) + 1.0);
      if (x < 0.25) {
        y = 1.0 - (x + x * y);
      } else {
        y = 0.5 - (x * y + (x - 0.5));
      }
    }
  } else if (absx < 1.25) {
    S = (absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) * ((absx -
      1.0) * ((absx - 1.0) * -0.0021663755948687908 + 0.035478304325618236) +
      -0.11089469428239668) + 0.31834661990116175) + -0.37220787603570132) +
                        0.41485611868374833) + -0.0023621185607526594;
    s = (absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) * ((absx -
      1.0) * ((absx - 1.0) * 0.011984499846799107 + 0.013637083912029051) +
      0.12617121980876164) + 0.071828654414196266) + 0.540397917702171) +
                        0.10642088040084423) + 1.0;
    if (x >= 0.0) {
      y = 0.15493708848953247 - S / s;
    } else {
      y = (S / s + 0.84506291151046753) + 1.0;
    }
  } else if (x < -6.0) {
    y = 2.0;
  } else if (x >= 28.0) {
    y = 0.0;
  } else {
    s = 1.0 / (absx * absx);
    if (absx < 2.8571414947509766) {
      R = s * (s * (s * (s * (s * (s * (s * -9.8143293441691455 +
        -81.2874355063066) + -184.60509290671104) + -162.39666946257347) +
                         -62.375332450326006) + -10.558626225323291) +
               -0.69385857270718176) + -0.0098649440348471482;
      S = s * (s * (s * (s * (s * (s * (s * (s * -0.0604244152148581 +
        6.5702497703192817) + 108.63500554177944) + 429.00814002756783) +
                  645.38727173326788) + 434.56587747522923) + 137.65775414351904)
               + 19.651271667439257) + 1.0;
    } else {
      R = s * (s * (s * (s * (s * (s * -483.5191916086514 + -1025.0951316110772)
                  + -637.56644336838963) + -160.63638485582192) +
                    -17.757954917754752) + -0.799283237680523) +
        -0.0098649429247001;
      S = s * (s * (s * (s * (s * (s * (s * -22.440952446585818 +
        474.52854120695537) + 2553.0504064331644) + 3199.8582195085955) +
                         1536.729586084437) + 325.79251299657392) +
               30.338060743482458) + 1.0;
    }

    if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
      s = frexp(absx, &eint);
    } else {
      s = absx;
      eint = 0;
    }

    s = std::floor(s * 2.097152E+6) / 2.097152E+6 * rt_powd_snf(2.0,
      static_cast<double>(eint));
    y = std::exp(-s * s - 0.5625) * std::exp((s - absx) * (s + absx) + R / S) /
      absx;
    if (x < 0.0) {
      y = 2.0 - y;
    }
  }

  return y;
}

//
// vp =   mvndnv( n, as, ch, bs, clg, ci, dci, xx, nv );
//   Transformed integrand for computation of MVN probabilities.
// Arguments    : double n
//                const emxArray_real_T *a
//                const emxArray_creal_T *ch
//                const emxArray_real_T *b
//                const emxArray_real_T *clg
//                double ci
//                double dci
//                const emxArray_real_T *x
//                double nv
//                emxArray_real_T *p
// Return Type  : void
//
static void mvndnv(double n, const emxArray_real_T *a, const emxArray_creal_T
                   *ch, const emxArray_real_T *b, const emxArray_real_T *clg,
                   double ci, double dci, const emxArray_real_T *x, double nv,
                   emxArray_real_T *p)
{
  emxArray_real_T *y;
  int i;
  int i1;
  int loop_ub;
  int b_loop_ub;
  emxArray_real_T *c;
  emxArray_real_T *dc;
  double li;
  double lf;
  emxArray_creal_T *s;
  emxArray_creal_T *tmp1;
  emxArray_real_T *bi;
  emxArray_real_T *r;
  emxArray_real_T *r1;
  emxArray_creal_T *b_y;
  int b_i;
  int i2;
  double z;
  double c_y;
  double mantissa;
  int nIterations;
  double b_x;
  boolean_T guard1 = false;
  int k;
  int exponent;
  double absx;
  int ch_re_tmp;
  double S;
  double R;
  int eint;
  emxInit_real_T(&y, 2);

  //  error estimate is 3 x standard error with ns samples.
  //
  //  end qscmvnv
  //
  i = static_cast<int>((n - 1.0));
  i1 = y->size[0] * y->size[1];
  y->size[0] = i;
  loop_ub = static_cast<int>(nv);
  y->size[1] = loop_ub;
  emxEnsureCapacity_real_T(y, i1);
  b_loop_ub = i * loop_ub;
  for (i = 0; i < b_loop_ub; i++) {
    y->data[i] = 0.0;
  }

  emxInit_real_T(&c, 2);
  i = c->size[0] * c->size[1];
  c->size[0] = 1;
  c->size[1] = loop_ub;
  emxEnsureCapacity_real_T(c, i);
  for (i = 0; i < loop_ub; i++) {
    c->data[i] = ci;
  }

  emxInit_real_T(&dc, 2);
  i = dc->size[0] * dc->size[1];
  dc->size[0] = 1;
  dc->size[1] = loop_ub;
  emxEnsureCapacity_real_T(dc, i);
  for (i = 0; i < loop_ub; i++) {
    dc->data[i] = dci;
  }

  i = p->size[0] * p->size[1];
  p->size[0] = 1;
  p->size[1] = dc->size[1];
  emxEnsureCapacity_real_T(p, i);
  b_loop_ub = dc->size[0] * dc->size[1];
  for (i = 0; i < b_loop_ub; i++) {
    p->data[i] = dc->data[i];
  }

  li = 2.0;
  lf = 1.0;
  i = static_cast<int>((n + -1.0));
  emxInit_creal_T(&s, 2);
  emxInit_creal_T(&tmp1, 2);
  emxInit_real_T(&bi, 2);
  emxInit_real_T(&r, 2);
  emxInit_real_T(&r1, 2);
  emxInit_creal_T(&b_y, 2);
  for (b_i = 0; b_i < i; b_i++) {
    i1 = dc->size[0] * dc->size[1];
    dc->size[0] = 1;
    b_loop_ub = c->size[1];
    dc->size[1] = c->size[1];
    emxEnsureCapacity_real_T(dc, i1);
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      dc->data[i1] = c->data[i1] + x->data[b_i + x->size[0] * i1] * dc->data[i1];
    }

    i1 = r->size[0] * r->size[1];
    r->size[0] = 1;
    i2 = dc->size[1];
    r->size[1] = dc->size[1];
    emxEnsureCapacity_real_T(r, i1);
    for (b_loop_ub = 0; b_loop_ub < i2; b_loop_ub++) {
      z = dc->data[b_loop_ub];
      if ((z >= 0.0) && (z <= 1.0)) {
        c_y = 2.0 * z;
        nIterations = 1;
        if (c_y == 0.0) {
          b_x = rtInf;
        } else if (c_y == 2.0) {
          b_x = rtMinusInf;
        } else {
          guard1 = false;
          if (c_y > 1.7) {
            z = std::sqrt(-std::log((2.0 - c_y) / 2.0));
            b_x = -(((1.641345311 * z + 3.429567803) * z + -1.624906493) * z +
                    -1.970840454) / ((1.6370678 * z + 3.5438892) * z + 1.0);
            guard1 = true;
          } else if (c_y < 0.3) {
            z = std::sqrt(0.69314718055994529 - std::log(c_y));
            if (c_y < 1.0947644252537633E-47) {
              if (c_y < 7.7532508072625747E-267) {
                nIterations = 3;
              } else {
                nIterations = 2;
              }
            }

            b_x = (((1.641345311 * z + 3.429567803) * z + -1.624906493) * z +
                   -1.970840454) / ((1.6370678 * z + 3.5438892) * z + 1.0);
            for (k = 0; k <= nIterations; k++) {
              z = -(eml_erfcore(b_x) - c_y) / (1.1283791670955126 * std::exp
                (-b_x * b_x));
              b_x -= z / (b_x * z + 1.0);
            }
          } else {
            z = (1.0 - c_y) * (1.0 - c_y);
            b_x = (1.0 - c_y) * (((-0.140543331 * z + 0.914624893) * z +
                                  -1.645349621) * z + 0.886226899) /
              ((((0.012229801 * z + -0.329097515) * z + 1.442710462) * z +
                -2.118377725) * z + 1.0);
            guard1 = true;
          }

          if (guard1) {
            if (c_y > 1.7) {
              for (k = 0; k < 2; k++) {
                z = (eml_erfcore(-b_x) - (2.0 - c_y)) / (1.1283791670955126 *
                  std::exp(-b_x * b_x));
                b_x -= z / (b_x * z + 1.0);
              }
            } else {
              for (k = 0; k < 2; k++) {
                // ========================== COPYRIGHT NOTICE ============================ 
                //  The algorithms for calculating ERF(X) and ERFC(X) are derived           
                //  from FDLIBM, which has the following notice:                            
                //                                                                          
                //  Copyright (C) 1993 by Sun Microsystems, Inc. All rights reserved.       
                //                                                                          
                //  Developed at SunSoft, a Sun Microsystems, Inc. business.                
                //  Permission to use, copy, modify, and distribute this                    
                //  software is freely granted, provided that this notice                   
                //  is preserved.                                                           
                // =============================    END    ================================ 
                absx = std::abs(b_x);
                if (rtIsNaN(b_x)) {
                  z = b_x;
                } else if (rtIsInf(b_x)) {
                  if (b_x < 0.0) {
                    z = -1.0;
                  } else {
                    z = 1.0;
                  }
                } else if (absx < 0.84375) {
                  if (absx < 3.7252902984619141E-9) {
                    if (absx < 2.8480945388892178E-306) {
                      z = 0.125 * (8.0 * b_x + 1.0270333367641007 * b_x);
                    } else {
                      z = b_x + 0.12837916709551259 * b_x;
                    }
                  } else {
                    z = b_x * b_x;
                    z = b_x + b_x * ((z * (z * (z * (z * -2.3763016656650163E-5
                      + -0.0057702702964894416) + -0.02848174957559851) +
                      -0.3250421072470015) + 0.12837916709551256) / (z * (z * (z
                      * (z * (z * -3.9602282787753681E-6 +
                              0.00013249473800432164) + 0.0050813062818757656) +
                      0.0650222499887673) + 0.39791722395915535) + 1.0));
                  }
                } else if (absx < 1.25) {
                  S = (absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) * ((absx -
                    1.0) * ((absx - 1.0) * ((absx - 1.0) *
                    -0.0021663755948687908 + 0.035478304325618236) +
                            -0.11089469428239668) + 0.31834661990116175) +
                    -0.37220787603570132) + 0.41485611868374833) +
                    -0.0023621185607526594;
                  z = (absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) * ((absx -
                    1.0) * ((absx - 1.0) * ((absx - 1.0) * 0.011984499846799107
                    + 0.013637083912029051) + 0.12617121980876164) +
                    0.071828654414196266) + 0.540397917702171) +
                                      0.10642088040084423) + 1.0;
                  if (b_x >= 0.0) {
                    z = S / z + 0.84506291151046753;
                  } else {
                    z = -0.84506291151046753 - S / z;
                  }
                } else if (absx > 6.0) {
                  if (b_x < 0.0) {
                    z = -1.0;
                  } else {
                    z = 1.0;
                  }
                } else {
                  z = 1.0 / (absx * absx);
                  if (absx < 2.8571434020996094) {
                    R = z * (z * (z * (z * (z * (z * (z * -9.8143293441691455 +
                      -81.2874355063066) + -184.60509290671104) +
                                -162.39666946257347) + -62.375332450326006) +
                                  -10.558626225323291) + -0.69385857270718176) +
                      -0.0098649440348471482;
                    S = z * (z * (z * (z * (z * (z * (z * (z *
                      -0.0604244152148581 + 6.5702497703192817) +
                      108.63500554177944) + 429.00814002756783) +
                                645.38727173326788) + 434.56587747522923) +
                                  137.65775414351904) + 19.651271667439257) +
                      1.0;
                  } else {
                    R = z * (z * (z * (z * (z * (z * -483.5191916086514 +
                      -1025.0951316110772) + -637.56644336838963) +
                                       -160.63638485582192) +
                                  -17.757954917754752) + -0.799283237680523) +
                      -0.0098649429247001;
                    S = z * (z * (z * (z * (z * (z * (z * -22.440952446585818 +
                      474.52854120695537) + 2553.0504064331644) +
                                3199.8582195085955) + 1536.729586084437) +
                                  325.79251299657392) + 30.338060743482458) +
                      1.0;
                  }

                  if (!rtIsNaN(absx)) {
                    mantissa = frexp(absx, &eint);
                    exponent = eint;
                  } else {
                    mantissa = absx;
                    exponent = 0;
                  }

                  z = std::floor(mantissa * 2.097152E+6) / 2.097152E+6 *
                    rt_powd_snf(2.0, static_cast<double>(exponent));
                  z = std::exp(-z * z - 0.5625) * std::exp((z - absx) * (z +
                    absx) + R / S) / absx;
                  if (b_x < 0.0) {
                    z--;
                  } else {
                    z = 1.0 - z;
                  }
                }

                z = (z - (1.0 - c_y)) / (1.1283791670955126 * std::exp(-b_x *
                  b_x));
                b_x -= z / (b_x * z + 1.0);
              }
            }
          }
        }

        r->data[b_loop_ub] = -1.4142135623730951 * b_x;
      } else {
        r->data[b_loop_ub] = rtNaN;
      }
    }

    b_loop_ub = r->size[1];
    for (i1 = 0; i1 < b_loop_ub; i1++) {
      y->data[b_i + y->size[0] * i1] = r->data[i1];
    }

    mantissa = clg->data[b_i + 1];
    lf += mantissa;
    if (lf < li) {
      i1 = c->size[0] * c->size[1];
      c->size[0] = 1;
      c->size[1] = 1;
      emxEnsureCapacity_real_T(c, i1);
      c->data[0] = 0.0;
      i1 = dc->size[0] * dc->size[1];
      dc->size[0] = 1;
      dc->size[1] = 1;
      emxEnsureCapacity_real_T(dc, i1);
      dc->data[0] = 1.0;
    } else {
      if (li > lf) {
        i1 = 0;
        i2 = 0;
      } else {
        i1 = static_cast<int>(li) - 1;
        i2 = static_cast<int>(lf);
      }

      b_loop_ub = y->size[1];
      nIterations = b_y->size[0] * b_y->size[1];
      b_y->size[0] = b_i + 1;
      b_y->size[1] = y->size[1];
      emxEnsureCapacity_creal_T(b_y, nIterations);
      for (nIterations = 0; nIterations < b_loop_ub; nIterations++) {
        for (k = 0; k <= b_i; k++) {
          b_y->data[k + b_y->size[0] * nIterations].re = y->data[k + y->size[0] *
            nIterations];
          b_y->data[k + b_y->size[0] * nIterations].im = 0.0;
        }
      }

      b_loop_ub = i2 - i1;
      i2 = s->size[0] * s->size[1];
      s->size[0] = b_loop_ub;
      s->size[1] = b_y->size[1];
      emxEnsureCapacity_creal_T(s, i2);
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        exponent = b_y->size[1];
        for (nIterations = 0; nIterations < exponent; nIterations++) {
          s->data[i2 + s->size[0] * nIterations].re = 0.0;
          s->data[i2 + s->size[0] * nIterations].im = 0.0;
          for (k = 0; k <= b_i; k++) {
            ch_re_tmp = i1 + i2;
            s->data[i2 + s->size[0] * nIterations].re += ch->data[ch_re_tmp +
              ch->size[0] * k].re * b_y->data[k + b_y->size[0] * nIterations].re
              - ch->data[ch_re_tmp + ch->size[0] * k].im * b_y->data[k +
              b_y->size[0] * nIterations].im;
            s->data[i2 + s->size[0] * nIterations].im += ch->data[ch_re_tmp +
              ch->size[0] * k].re * b_y->data[k + b_y->size[0] * nIterations].im
              + ch->data[ch_re_tmp + ch->size[0] * k].im * b_y->data[k +
              b_y->size[0] * nIterations].re;
          }
        }
      }

      //  original code: ai =          max( max( a(li:lf)*on - s, [], 1 ), -36.0 );  
      if (li > lf) {
        i1 = 0;
        i2 = 0;
      } else {
        i1 = static_cast<int>(li) - 1;
        i2 = static_cast<int>(lf);
      }

      b_loop_ub = i2 - i1;
      i2 = tmp1->size[0] * tmp1->size[1];
      tmp1->size[0] = b_loop_ub;
      tmp1->size[1] = loop_ub;
      emxEnsureCapacity_creal_T(tmp1, i2);
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        for (nIterations = 0; nIterations < loop_ub; nIterations++) {
          tmp1->data[i2 + tmp1->size[0] * nIterations].re = a->data[i1 + i2] -
            s->data[i2 + s->size[0] * nIterations].re;
          tmp1->data[i2 + tmp1->size[0] * nIterations].im = 0.0 - s->data[i2 +
            s->size[0] * nIterations].im;
        }
      }

      i1 = dc->size[0] * dc->size[1];
      dc->size[0] = 1;
      b_loop_ub = tmp1->size[1];
      dc->size[1] = tmp1->size[1];
      emxEnsureCapacity_real_T(dc, i1);
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        dc->data[i1] = 0.0;
      }

      i1 = tmp1->size[1];
      for (b_loop_ub = 0; b_loop_ub < i1; b_loop_ub++) {
        S = -36.0;
        i2 = tmp1->size[0];
        for (nIterations = 0; nIterations < i2; nIterations++) {
          if (tmp1->data[nIterations + tmp1->size[0] * b_loop_ub].re > S) {
            S = tmp1->data[nIterations + tmp1->size[0] * b_loop_ub].re;
          }
        }

        dc->data[b_loop_ub] = S;
      }

      //  bi = max( ai, min( min( b(li:lf)*on - s, [], 1 ),   9 ) );
      if (li > lf) {
        i1 = 0;
        i2 = 0;
      } else {
        i1 = static_cast<int>(li) - 1;
        i2 = static_cast<int>(lf);
      }

      b_loop_ub = i2 - i1;
      i2 = tmp1->size[0] * tmp1->size[1];
      tmp1->size[0] = b_loop_ub;
      tmp1->size[1] = loop_ub;
      emxEnsureCapacity_creal_T(tmp1, i2);
      for (i2 = 0; i2 < b_loop_ub; i2++) {
        for (nIterations = 0; nIterations < loop_ub; nIterations++) {
          tmp1->data[i2 + tmp1->size[0] * nIterations].re = b->data[i1 + i2] -
            s->data[i2 + s->size[0] * nIterations].re;
          tmp1->data[i2 + tmp1->size[0] * nIterations].im = 0.0 - s->data[i2 +
            s->size[0] * nIterations].im;
        }
      }

      i1 = bi->size[0] * bi->size[1];
      bi->size[0] = 1;
      b_loop_ub = tmp1->size[1];
      bi->size[1] = tmp1->size[1];
      emxEnsureCapacity_real_T(bi, i1);
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        bi->data[i1] = 0.0;
      }

      i1 = tmp1->size[1];
      for (b_loop_ub = 0; b_loop_ub < i1; b_loop_ub++) {
        S = 9.0;
        i2 = tmp1->size[0];
        for (nIterations = 0; nIterations < i2; nIterations++) {
          if (tmp1->data[nIterations + tmp1->size[0] * b_loop_ub].re < S) {
            S = tmp1->data[nIterations + tmp1->size[0] * b_loop_ub].re;
          }
        }

        z = dc->data[b_loop_ub];
        if (z > S) {
          S = z;
        }

        bi->data[b_loop_ub] = S;
      }

      i1 = c->size[0] * c->size[1];
      c->size[0] = 1;
      i2 = dc->size[1];
      c->size[1] = dc->size[1];
      emxEnsureCapacity_real_T(c, i1);
      for (b_loop_ub = 0; b_loop_ub < i2; b_loop_ub++) {
        c_y = eml_erfcore(-dc->data[b_loop_ub] / 1.4142135623730951);
        c->data[b_loop_ub] = 0.5 * c_y;
      }

      i1 = r1->size[0] * r1->size[1];
      r1->size[0] = 1;
      i2 = bi->size[1];
      r1->size[1] = bi->size[1];
      emxEnsureCapacity_real_T(r1, i1);
      for (b_loop_ub = 0; b_loop_ub < i2; b_loop_ub++) {
        c_y = eml_erfcore(-bi->data[b_loop_ub] / 1.4142135623730951);
        r1->data[b_loop_ub] = 0.5 * c_y;
      }

      i1 = dc->size[0] * dc->size[1];
      dc->size[0] = 1;
      dc->size[1] = r1->size[1];
      emxEnsureCapacity_real_T(dc, i1);
      b_loop_ub = r1->size[0] * r1->size[1];
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        dc->data[i1] = r1->data[i1] - c->data[i1];
      }

      i1 = p->size[0] * p->size[1];
      i2 = p->size[0] * p->size[1];
      p->size[0] = 1;
      emxEnsureCapacity_real_T(p, i2);
      b_loop_ub = i1 - 1;
      for (i1 = 0; i1 <= b_loop_ub; i1++) {
        p->data[i1] *= dc->data[i1];
      }
    }

    li += mantissa;
  }

  emxFree_creal_T(&b_y);
  emxFree_real_T(&r1);
  emxFree_real_T(&r);
  emxFree_real_T(&bi);
  emxFree_creal_T(&tmp1);
  emxFree_creal_T(&s);
  emxFree_real_T(&dc);
  emxFree_real_T(&c);
  emxFree_real_T(&y);
}

//
// Arguments    : double u0
//                double u1
// Return Type  : double
//
static double rt_powd_snf(double u0, double u1)
{
  double y;
  double d;
  double d1;
  if (rtIsNaN(u0) || rtIsNaN(u1)) {
    y = rtNaN;
  } else {
    d = std::abs(u0);
    d1 = std::abs(u1);
    if (rtIsInf(u1)) {
      if (d == 1.0) {
        y = 1.0;
      } else if (d > 1.0) {
        if (u1 > 0.0) {
          y = rtInf;
        } else {
          y = 0.0;
        }
      } else if (u1 > 0.0) {
        y = 0.0;
      } else {
        y = rtInf;
      }
    } else if (d1 == 0.0) {
      y = 1.0;
    } else if (d1 == 1.0) {
      if (u1 > 0.0) {
        y = u0;
      } else {
        y = 1.0 / u0;
      }
    } else if (u1 == 2.0) {
      y = u0 * u0;
    } else if ((u1 == 0.5) && (u0 >= 0.0)) {
      y = std::sqrt(u0);
    } else if ((u0 < 0.0) && (u1 > std::floor(u1))) {
      y = rtNaN;
    } else {
      y = pow(u0, u1);
    }
  }

  return y;
}

//
// [ P E ] = QSCMVNV( M, R, A, CN, B )
//     uses a randomized quasi-random rule with m points to estimate an
//     MVN probability for positive semi-definite covariance matrix r,
//     with constraints a < cn*x < b. If r is nxn and cn is kxn, then
//     a and b must be column k-vectors.
//    Probability p is output with error estimate e.
//     Example use:
//      r = [ 4 3 2 1; 3 5 -1 1; 2 -1 4 2; 1 1 2 5 ];
//      a = [ -inf 1 -5 ]'; b = [ 3 inf 4 ]';
//      cn = [ 1 2 3 -2; 2 4 1 2; -2 3 4 1 ];
//      [ p e ] = qscmvnv( 5000, r, a, cn, b ); disp([ p e ])
//
//   This function uses an algorithm given in the paper by Alan Genz:
//    "Numerical Computation of Multivariate Normal Probabilities", in
//      J. of Computational and Graphical Stat., 1(1992), 141-149.
//   The primary references for the numerical integration are
//    "On a Number-Theoretical Integration Method"
//      H. Niederreiter, Aequationes Mathematicae, 8(1972), 304-11, and
//    "Randomization of Number Theoretic Methods for Multiple Integration"
//      R. Cranley and T.N.L. Patterson, SIAM J Numer Anal, 13(1976), 904-14.
//
//    Alan Genz is the author of this function and following Matlab functions.
//           Alan Genz, WSU Math, PO Box 643113, Pullman, WA 99164-3113
//           Email : AlanGenz@wsu.edu
//
//
//  Copyright (C) 2014, Alan Genz,  All rights reserved.
//
//  Redistribution and use in source and binary forms, with or without
//  modification, are permitted provided the following conditions are met:
//    1. Redistributions of source code must retain the above copyright
//       notice, this list of conditions and the following disclaimer.
//    2. Redistributions in binary form must reproduce the above copyright
//       notice, this list of conditions and the following disclaimer in
//       the documentation and/or other materials provided with the
//       distribution.
//    3. The contributor name(s) may not be used to endorse or promote
//       products derived from this software without specific prior
//       written permission.
//  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
//  "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
//  LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS
//  FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
//  COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
//  INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
//  BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
//  OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND
//  ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR
//  TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF USE
//  OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
//  Initialization
// Arguments    : double m
//                const emxArray_real_T *r
//                const emxArray_real_T *a
//                const emxArray_real_T *cn
//                const emxArray_real_T *b
//                double *p
//                double *e
// Return Type  : void
//
void qscmvnvR(double m, const emxArray_real_T *r, const emxArray_real_T *a,
              const emxArray_real_T *cn, const emxArray_real_T *b, double *p,
              double *e)
{
  unsigned int b_r;
  int nx;
  emxArray_real_T *as;
  int i;
  int loop_ub;
  emxArray_real_T *bs;
  emxArray_creal_T *ch;
  emxArray_real_T *clg;
  double n;
  double d;
  double ci;
  double dci;
  double nv;
  emxArray_real_T *ps;
  int k;
  unsigned int np;
  emxArray_real_T *xx;
  emxArray_boolean_T *isp;
  int i1;
  unsigned int kmax;
  unsigned int b_k;
  int b_loop_ub;
  unsigned int nleft;
  unsigned int j;
  emxArray_real_T *b_ps;
  emxArray_real_T *b_a;
  emxArray_real_T *tmpDo;
  emxArray_real_T *c_a;
  emxArray_real_T *b_b;
  emxArray_real_T *x;
  int b_i;
  double b_d;
  if (isInitialized_qscmvnvR == false) {
    qscmvnvR_initialize();
  }

  std::memset(&state[0], 0, 625U * sizeof(unsigned int));
  b_r = 5489U;
  state[0] = 5489U;
  for (nx = 0; nx < 623; nx++) {
    b_r = ((b_r ^ b_r >> 30U) * 1812433253U + nx) + 1U;
    state[nx + 1] = b_r;
  }

  emxInit_real_T(&as, 1);
  state[624] = 624U;

  // as = zeros(size(a,1),1);
  i = as->size[0];
  as->size[0] = a->size[0];
  emxEnsureCapacity_real_T(as, i);
  loop_ub = a->size[0];
  for (i = 0; i < loop_ub; i++) {
    as->data[i] = a->data[i];
  }

  emxInit_real_T(&bs, 1);
  i = bs->size[0];
  bs->size[0] = b->size[0];
  emxEnsureCapacity_real_T(bs, i);
  loop_ub = b->size[0];
  for (i = 0; i < loop_ub; i++) {
    bs->data[i] = b->data[i];
  }

  emxInit_creal_T(&ch, 2);
  emxInit_real_T(&clg, 2);
  chlsrt(r, as, cn, bs, ch, clg, &n);
  d = eml_erfcore(-as->data[0] / 1.4142135623730951);
  ci = 0.5 * d;
  d = eml_erfcore(-bs->data[0] / 1.4142135623730951);
  dci = 0.5 * d - ci;
  *p = 0.0;
  *e = 0.0;
  nv = m / 12.0;
  if ((nv < 1.0) || rtIsNaN(nv)) {
    nv = 1.0;
  }

  emxInit_real_T(&ps, 2);
  nv = std::floor(nv);
  d = 5.0 * (n + 1.0) * std::log(n + 1.0) / 4.0;
  if (d < 2.0) {
    ps->size[0] = 1;
    ps->size[1] = 0;
  } else if (d < 3.0) {
    i = ps->size[0] * ps->size[1];
    ps->size[0] = 1;
    ps->size[1] = 1;
    emxEnsureCapacity_real_T(ps, i);
    ps->data[0] = 2.0;
  } else {
    d = std::floor(d);
    if (d < 4.294967296E+9) {
      b_r = static_cast<unsigned int>(d);
    } else if (d >= 4.294967296E+9) {
      b_r = MAX_uint32_T;
    } else {
      b_r = 0U;
    }

    np = b_r >> 1;
    if (b_r - (np << 1) > 0U) {
      np++;
    }

    emxInit_boolean_T(&isp, 2);
    i = isp->size[0] * isp->size[1];
    isp->size[0] = 1;
    isp->size[1] = static_cast<int>(np);
    emxEnsureCapacity_boolean_T(isp, i);
    loop_ub = static_cast<int>(np);
    for (i = 0; i < loop_ub; i++) {
      isp->data[i] = true;
    }

    kmax = intsqrt(b_r);
    b_k = 3U;
    nleft = np;
    while (b_k <= kmax) {
      if (isp->data[static_cast<int>(((b_k + 1U) >> 1)) - 1]) {
        b_r = (b_k * b_k + 1U) >> 1;
        for (j = b_r; j <= np; j += b_k) {
          i = static_cast<int>(j) - 1;
          if (isp->data[i]) {
            isp->data[i] = false;
            nleft--;
          }
        }
      }

      b_k += 2U;
    }

    i = ps->size[0] * ps->size[1];
    ps->size[0] = 1;
    ps->size[1] = static_cast<int>(nleft);
    emxEnsureCapacity_real_T(ps, i);
    b_r = 1U;
    ps->data[0] = 2.0;
    for (b_k = 2U; b_k <= np; b_k++) {
      if (isp->data[static_cast<int>(b_k) - 1] && (b_r < 203280221U)) {
        b_r++;
        ps->data[static_cast<int>(b_r) - 1] = (b_k << 1) - 1U;
      }
    }

    emxFree_boolean_T(&isp);
  }

  nx = ps->size[1];
  for (k = 0; k < nx; k++) {
    ps->data[k] = std::sqrt(ps->data[k]);
  }

  emxInit_real_T(&xx, 2);
  if (1.0 > n - 1.0) {
    loop_ub = 0;
  } else {
    loop_ub = static_cast<int>((n - 1.0));
  }

  //  Richtmyer generators
  //
  //  Randomization loop for ns samples
  //
  i = static_cast<int>((n - 1.0));
  i1 = xx->size[0] * xx->size[1];
  xx->size[0] = i;
  b_loop_ub = static_cast<int>(nv);
  xx->size[1] = b_loop_ub;
  emxEnsureCapacity_real_T(xx, i1);
  nx = i * b_loop_ub;
  for (i = 0; i < nx; i++) {
    xx->data[i] = 0.0;
  }

  emxInit_real_T(&b_ps, 2);
  emxInit_real_T(&b_a, 2);

  //  added
  emxInit_real_T(&tmpDo, 2);
  emxInit_real_T(&c_a, 1);
  emxInit_real_T(&b_b, 2);
  emxInit_real_T(&x, 2);
  for (b_i = 0; b_i < 12; b_i++) {
    if (rtIsInf(nv) && (1.0 == nv)) {
      i = tmpDo->size[0] * tmpDo->size[1];
      tmpDo->size[0] = 1;
      tmpDo->size[1] = 1;
      emxEnsureCapacity_real_T(tmpDo, i);
      tmpDo->data[0] = rtNaN;
    } else {
      i = tmpDo->size[0] * tmpDo->size[1];
      tmpDo->size[0] = 1;
      nx = static_cast<int>((nv - 1.0));
      tmpDo->size[1] = nx + 1;
      emxEnsureCapacity_real_T(tmpDo, i);
      for (i = 0; i <= nx; i++) {
        tmpDo->data[i] = static_cast<double>(i) + 1.0;
      }
    }

    b_rand(n - 1.0, c_a);
    i = b_ps->size[0] * b_ps->size[1];
    b_ps->size[0] = loop_ub;
    b_ps->size[1] = tmpDo->size[1];
    emxEnsureCapacity_real_T(b_ps, i);
    for (i = 0; i < loop_ub; i++) {
      nx = tmpDo->size[1];
      for (i1 = 0; i1 < nx; i1++) {
        b_ps->data[i + b_ps->size[0] * i1] = ps->data[i] * tmpDo->data[i1];
      }
    }

    i = b_a->size[0] * b_a->size[1];
    b_a->size[0] = c_a->size[0];
    b_a->size[1] = b_loop_ub;
    emxEnsureCapacity_real_T(b_a, i);
    nx = c_a->size[0];
    for (i = 0; i < nx; i++) {
      for (i1 = 0; i1 < b_loop_ub; i1++) {
        b_a->data[i + b_a->size[0] * i1] = c_a->data[i];
      }
    }

    i = x->size[0] * x->size[1];
    x->size[0] = b_ps->size[0];
    x->size[1] = b_ps->size[1];
    emxEnsureCapacity_real_T(x, i);
    nx = b_ps->size[0] * b_ps->size[1];
    for (i = 0; i < nx; i++) {
      x->data[i] = b_ps->data[i] + b_a->data[i];
    }

    i = x->size[0];
    i1 = b_b->size[0] * b_b->size[1];
    b_b->size[0] = i;
    nx = x->size[1];
    b_b->size[1] = nx;
    emxEnsureCapacity_real_T(b_b, i1);
    nx *= i;
    for (k = 0; k < nx; k++) {
      if (rtIsNaN(x->data[k]) || rtIsInf(x->data[k])) {
        d = rtNaN;
      } else if (x->data[k] == 0.0) {
        d = 0.0;
      } else {
        d = std::fmod(x->data[k], 1.0);
        if (d == 0.0) {
          d = 0.0;
        } else {
          if (x->data[k] < 0.0) {
            d++;
          }
        }
      }

      b_b->data[k] = d;
    }

    i = x->size[0] * x->size[1];
    x->size[0] = b_b->size[0];
    x->size[1] = b_b->size[1];
    emxEnsureCapacity_real_T(x, i);
    nx = b_b->size[0] * b_b->size[1];
    for (i = 0; i < nx; i++) {
      x->data[i] = 2.0 * b_b->data[i] - 1.0;
    }

    nx = x->size[0] * x->size[1];
    i = b_ps->size[0] * b_ps->size[1];
    b_ps->size[0] = x->size[0];
    b_ps->size[1] = x->size[1];
    emxEnsureCapacity_real_T(b_ps, i);
    for (k = 0; k < nx; k++) {
      b_ps->data[k] = std::abs(x->data[k]);
    }

    nx = b_ps->size[1];
    for (i = 0; i < nx; i++) {
      k = b_ps->size[0];
      for (i1 = 0; i1 < k; i1++) {
        xx->data[i1 + xx->size[0] * i] = b_ps->data[i1 + b_ps->size[0] * i];
      }
    }

    mvndnv(n, as, ch, bs, clg, ci, dci, xx, nv, tmpDo);
    i = tmpDo->size[1];
    if (tmpDo->size[1] == 0) {
      d = 0.0;
    } else {
      d = tmpDo->data[0];
      for (k = 2; k <= i; k++) {
        d += tmpDo->data[k - 1];
      }
    }

    d = (d / static_cast<double>(tmpDo->size[1]) - *p) / (static_cast<double>
      (b_i) + 1.0);
    *p += d;
    b_d = std::abs(d);
    if (b_d > 0.0) {
      d = *e / d;
      *e = b_d * std::sqrt(d * d * ((static_cast<double>(b_i) + 1.0) - 2.0) / (
        static_cast<double>(b_i) + 1.0) + 1.0);
    } else {
      if (b_i + 1 > 1) {
        *e *= std::sqrt(((static_cast<double>(b_i) + 1.0) - 2.0) / (static_cast<
          double>(b_i) + 1.0));
      }
    }
  }

  emxFree_real_T(&b_a);
  emxFree_real_T(&b_ps);
  emxFree_real_T(&x);
  emxFree_real_T(&b_b);
  emxFree_real_T(&c_a);
  emxFree_real_T(&clg);
  emxFree_real_T(&bs);
  emxFree_creal_T(&ch);
  emxFree_real_T(&as);
  emxFree_real_T(&tmpDo);
  emxFree_real_T(&xx);
  emxFree_real_T(&ps);
  *e *= 3.0;
}

//
// File trailer for qscmvnvR.cpp
//
// [EOF]
//
