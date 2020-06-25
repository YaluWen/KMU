//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: cnmatrix.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 17:21:40
//

// Include Files
#include "cnmatrix.h"
#include "cnmatrix_data.h"
#include "cnmatrix_emxutil.h"
#include "cnmatrix_initialize.h"
#include "cnmatrix_rtwutil.h"
#include "erfcinv.h"
#include "rt_nonfinite.h"
#include "sortLE.h"
#include <cmath>

// Function Definitions

//
// NoKernel=3;
// pobs=0.05;
// CovVar=eye(NoKernel);
// CovVar(1,1)=1;
// CovVar(1,2)=0.4;
// CovVar(2,1)=0.4;
// Arguments    : double NoKernel
//                double pobs
//                const emxArray_real_T *V
//                const emxArray_real_T *D
//                emxArray_real_T *r1
//                emxArray_real_T *a1
//                emxArray_real_T *cn1
//                emxArray_real_T *b1
// Return Type  : void
//
void cnmatrix(double NoKernel, double pobs, const emxArray_real_T *V, const
              emxArray_real_T *D, emxArray_real_T *r1, emxArray_real_T *a1,
              emxArray_real_T *cn1, emxArray_real_T *b1)
{
  double t;
  int m;
  int i;
  int qEnd;
  int k;
  emxArray_real_T *cn;
  int i1;
  int vstride;
  emxArray_int8_T *r;
  emxArray_int8_T *ycol;
  int b_i;
  int iy;
  emxArray_boolean_T *b_r1;
  unsigned int outsize_idx_0;
  int i2;
  int ibcol;
  int j;
  emxArray_int32_T *idx;
  boolean_T exitg1;
  emxArray_real_T *b_cn;
  emxArray_real_T *y;
  emxArray_int32_T *col;
  int n;
  emxArray_int32_T *iwork;
  emxArray_real_T *b_y;
  if (isInitialized_cnmatrix == false) {
    cnmatrix_initialize();
  }

  if (NoKernel < 0.0) {
    t = 0.0;
  } else {
    t = NoKernel;
  }

  m = static_cast<int>(t);
  i = r1->size[0] * r1->size[1];
  r1->size[0] = static_cast<int>(t);
  r1->size[1] = static_cast<int>(t);
  emxEnsureCapacity_real_T(r1, i);
  qEnd = static_cast<int>(t) * static_cast<int>(t);
  for (i = 0; i < qEnd; i++) {
    r1->data[i] = 0.0;
  }

  if (static_cast<int>(t) > 0) {
    for (k = 0; k < m; k++) {
      r1->data[k + r1->size[0] * k] = 1.0;
    }
  }

  emxInit_real_T(&cn, 2);
  t = rt_powd_snf(3.0, NoKernel);
  i = static_cast<int>(t);
  i1 = cn->size[0] * cn->size[1];
  cn->size[0] = i;
  vstride = static_cast<int>(NoKernel);
  cn->size[1] = vstride;
  emxEnsureCapacity_real_T(cn, i1);
  qEnd = i * vstride;
  for (i = 0; i < qEnd; i++) {
    cn->data[i] = 0.0;
  }

  emxInit_int8_T(&r, 1);
  emxInit_int8_T(&ycol, 1);
  for (b_i = 0; b_i < vstride; b_i++) {
    iy = static_cast<int>(rt_powd_snf(3.0, (static_cast<double>(b_i) + 1.0) -
      1.0));
    qEnd = static_cast<int>(rt_powd_snf(3.0, (static_cast<double>(b_i) + 1.0) -
      1.0));
    i = ycol->size[0];
    ycol->size[0] = (iy + iy) + qEnd;
    emxEnsureCapacity_int8_T(ycol, i);
    for (i = 0; i < iy; i++) {
      ycol->data[i] = 1;
    }

    for (i = 0; i < iy; i++) {
      ycol->data[i + iy] = 0;
    }

    for (i = 0; i < qEnd; i++) {
      ycol->data[(i + iy) + iy] = -1;
    }

    i = ycol->size[0];
    i1 = static_cast<int>((t / rt_powd_snf(3.0, static_cast<double>(b_i) + 1.0)));
    iy = r->size[0];
    r->size[0] = ycol->size[0] * i1;
    emxEnsureCapacity_int8_T(r, iy);
    for (iy = 0; iy < i1; iy++) {
      ibcol = iy * i;
      for (k = 0; k < i; k++) {
        r->data[ibcol + k] = ycol->data[k];
      }
    }

    qEnd = r->size[0];
    for (i = 0; i < qEnd; i++) {
      cn->data[i + cn->size[0] * b_i] = r->data[i];
    }
  }

  emxFree_int8_T(&r);
  emxInit_boolean_T(&b_r1, 1);

  //  Reduce redundency %
  //  remove zero rows %
  outsize_idx_0 = static_cast<unsigned int>(cn->size[0]);
  i = b_r1->size[0];
  b_r1->size[0] = static_cast<int>(outsize_idx_0);
  emxEnsureCapacity_boolean_T(b_r1, i);
  qEnd = static_cast<int>(outsize_idx_0);
  for (i = 0; i < qEnd; i++) {
    b_r1->data[i] = false;
  }

  vstride = cn->size[0];
  i2 = (cn->size[1] - 1) * cn->size[0];
  iy = -1;
  ibcol = 0;
  for (j = 0; j < vstride; j++) {
    ibcol++;
    i2++;
    iy++;
    qEnd = ibcol;
    exitg1 = false;
    while ((!exitg1) && ((vstride > 0) && (qEnd <= i2))) {
      if (static_cast<signed char>(cn->data[qEnd - 1]) == 0) {
        qEnd += vstride;
      } else {
        b_r1->data[iy] = true;
        exitg1 = true;
      }
    }
  }

  ibcol = b_r1->size[0] - 1;
  iy = 0;
  for (b_i = 0; b_i <= ibcol; b_i++) {
    if (b_r1->data[b_i]) {
      iy++;
    }
  }

  emxInit_int32_T(&idx, 1);
  i = idx->size[0];
  idx->size[0] = iy;
  emxEnsureCapacity_int32_T(idx, i);
  iy = 0;
  for (b_i = 0; b_i <= ibcol; b_i++) {
    if (b_r1->data[b_i]) {
      idx->data[iy] = b_i + 1;
      iy++;
    }
  }

  emxFree_boolean_T(&b_r1);
  emxInit_real_T(&b_cn, 2);
  qEnd = cn->size[1] - 1;
  i = b_cn->size[0] * b_cn->size[1];
  b_cn->size[0] = idx->size[0];
  b_cn->size[1] = cn->size[1];
  emxEnsureCapacity_real_T(b_cn, i);
  for (i = 0; i <= qEnd; i++) {
    iy = idx->size[0];
    for (i1 = 0; i1 < iy; i1++) {
      b_cn->data[i1 + b_cn->size[0] * i] = cn->data[(idx->data[i1] + cn->size[0]
        * i) - 1];
    }
  }

  i = cn->size[0] * cn->size[1];
  cn->size[0] = b_cn->size[0];
  cn->size[1] = b_cn->size[1];
  emxEnsureCapacity_real_T(cn, i);
  qEnd = b_cn->size[0] * b_cn->size[1];
  for (i = 0; i < qEnd; i++) {
    cn->data[i] = b_cn->data[i];
  }

  emxInit_real_T(&y, 2);

  //  keep only one with A=-A %
  if (cn->size[1] < 1) {
    y->size[0] = 1;
    y->size[1] = 0;
  } else {
    i = y->size[0] * y->size[1];
    y->size[0] = 1;
    y->size[1] = static_cast<int>((static_cast<double>(cn->size[1]) - 1.0)) + 1;
    emxEnsureCapacity_real_T(y, i);
    qEnd = static_cast<int>((static_cast<double>(cn->size[1]) - 1.0));
    for (i = 0; i <= qEnd; i++) {
      y->data[i] = static_cast<double>(i) + 1.0;
    }
  }

  emxInit_int32_T(&col, 2);
  i = col->size[0] * col->size[1];
  col->size[0] = 1;
  col->size[1] = y->size[1];
  emxEnsureCapacity_int32_T(col, i);
  qEnd = y->size[0] * y->size[1];
  for (i = 0; i < qEnd; i++) {
    col->data[i] = static_cast<int>(y->data[i]);
  }

  emxFree_real_T(&y);
  n = cn->size[0] + 1;
  i = idx->size[0];
  idx->size[0] = cn->size[0];
  emxEnsureCapacity_int32_T(idx, i);
  qEnd = cn->size[0];
  for (i = 0; i < qEnd; i++) {
    idx->data[i] = 0;
  }

  if ((cn->size[0] == 0) || (cn->size[1] == 0)) {
    for (k = 0; k <= n - 2; k++) {
      idx->data[k] = k + 1;
    }
  } else {
    emxInit_int32_T(&iwork, 1);
    i = iwork->size[0];
    iwork->size[0] = cn->size[0];
    emxEnsureCapacity_int32_T(iwork, i);
    i = cn->size[0] - 1;
    for (k = 1; k <= i; k += 2) {
      if (sortLE(cn, col, k, k + 1)) {
        idx->data[k - 1] = k;
        idx->data[k] = k + 1;
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    if ((cn->size[0] & 1) != 0) {
      idx->data[cn->size[0] - 1] = cn->size[0];
    }

    b_i = 2;
    while (b_i < n - 1) {
      i2 = b_i << 1;
      j = 1;
      for (iy = b_i + 1; iy < n; iy = qEnd + b_i) {
        ibcol = j;
        vstride = iy;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }

        k = 0;
        m = qEnd - j;
        while (k + 1 <= m) {
          i = idx->data[vstride - 1];
          i1 = idx->data[ibcol - 1];
          if (sortLE(cn, col, i1, i)) {
            iwork->data[k] = i1;
            ibcol++;
            if (ibcol == iy) {
              while (vstride < qEnd) {
                k++;
                iwork->data[k] = idx->data[vstride - 1];
                vstride++;
              }
            }
          } else {
            iwork->data[k] = i;
            vstride++;
            if (vstride == qEnd) {
              while (ibcol < iy) {
                k++;
                iwork->data[k] = idx->data[ibcol - 1];
                ibcol++;
              }
            }
          }

          k++;
        }

        for (k = 0; k < m; k++) {
          idx->data[(j + k) - 1] = iwork->data[k];
        }

        j = qEnd;
      }

      b_i = i2;
    }

    emxFree_int32_T(&iwork);
  }

  emxFree_int32_T(&col);
  m = cn->size[0] - 1;
  n = cn->size[1];
  i = ycol->size[0];
  ycol->size[0] = cn->size[0];
  emxEnsureCapacity_int8_T(ycol, i);
  for (j = 0; j < n; j++) {
    for (b_i = 0; b_i <= m; b_i++) {
      ycol->data[b_i] = static_cast<signed char>(cn->data[(idx->data[b_i] +
        cn->size[0] * j) - 1]);
    }

    for (b_i = 0; b_i <= m; b_i++) {
      cn->data[b_i + cn->size[0] * j] = ycol->data[b_i];
    }
  }

  emxFree_int8_T(&ycol);
  emxFree_int32_T(&idx);
  t = static_cast<double>(cn->size[0]) / 2.0;
  if (1.0 > t) {
    qEnd = 0;
  } else {
    qEnd = static_cast<int>(t);
  }

  iy = cn->size[1] - 1;
  i = b_cn->size[0] * b_cn->size[1];
  b_cn->size[0] = qEnd;
  b_cn->size[1] = cn->size[1];
  emxEnsureCapacity_real_T(b_cn, i);
  for (i = 0; i <= iy; i++) {
    for (i1 = 0; i1 < qEnd; i1++) {
      b_cn->data[i1 + b_cn->size[0] * i] = cn->data[i1 + cn->size[0] * i];
    }
  }

  i = cn->size[0] * cn->size[1];
  cn->size[0] = b_cn->size[0];
  cn->size[1] = b_cn->size[1];
  emxEnsureCapacity_real_T(cn, i);
  qEnd = b_cn->size[0] * b_cn->size[1];
  for (i = 0; i < qEnd; i++) {
    cn->data[i] = b_cn->data[i];
  }

  emxFree_real_T(&b_cn);
  emxInit_real_T(&b_y, 2);

  // [V,D] = eig(CovVar);
  if ((cn->size[1] == 1) || (V->size[0] == 1)) {
    i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = cn->size[0];
    b_y->size[1] = V->size[1];
    emxEnsureCapacity_real_T(b_y, i);
    qEnd = cn->size[0];
    for (i = 0; i < qEnd; i++) {
      iy = V->size[1];
      for (i1 = 0; i1 < iy; i1++) {
        b_y->data[i + b_y->size[0] * i1] = 0.0;
        ibcol = cn->size[1];
        for (vstride = 0; vstride < ibcol; vstride++) {
          b_y->data[i + b_y->size[0] * i1] += cn->data[i + cn->size[0] * vstride]
            * V->data[vstride + V->size[0] * i1];
        }
      }
    }
  } else {
    m = cn->size[0];
    iy = cn->size[1];
    n = V->size[1];
    i = b_y->size[0] * b_y->size[1];
    b_y->size[0] = cn->size[0];
    b_y->size[1] = V->size[1];
    emxEnsureCapacity_real_T(b_y, i);
    for (j = 0; j < n; j++) {
      ibcol = j * m;
      qEnd = j * iy;
      for (b_i = 0; b_i < m; b_i++) {
        b_y->data[ibcol + b_i] = 0.0;
      }

      for (k = 0; k < iy; k++) {
        vstride = k * m;
        t = V->data[qEnd + k];
        for (b_i = 0; b_i < m; b_i++) {
          i = ibcol + b_i;
          b_y->data[i] += t * static_cast<double>(static_cast<signed char>
            (cn->data[vstride + b_i]));
        }
      }
    }
  }

  i = cn->size[0] * cn->size[1];
  cn->size[0] = D->size[0];
  cn->size[1] = D->size[1];
  emxEnsureCapacity_real_T(cn, i);
  qEnd = D->size[0] * D->size[1];
  for (i = 0; i < qEnd; i++) {
    cn->data[i] = D->data[i];
  }

  iy = D->size[0] * D->size[1];
  for (k = 0; k < iy; k++) {
    cn->data[k] = std::sqrt(cn->data[k]);
  }

  if ((b_y->size[1] == 1) || (cn->size[0] == 1)) {
    i = cn1->size[0] * cn1->size[1];
    cn1->size[0] = b_y->size[0];
    cn1->size[1] = cn->size[1];
    emxEnsureCapacity_real_T(cn1, i);
    qEnd = b_y->size[0];
    for (i = 0; i < qEnd; i++) {
      iy = cn->size[1];
      for (i1 = 0; i1 < iy; i1++) {
        cn1->data[i + cn1->size[0] * i1] = 0.0;
        ibcol = b_y->size[1];
        for (vstride = 0; vstride < ibcol; vstride++) {
          cn1->data[i + cn1->size[0] * i1] += b_y->data[i + b_y->size[0] *
            vstride] * cn->data[vstride + cn->size[0] * i1];
        }
      }
    }
  } else {
    m = b_y->size[0];
    iy = b_y->size[1];
    n = cn->size[1];
    i = cn1->size[0] * cn1->size[1];
    cn1->size[0] = b_y->size[0];
    cn1->size[1] = cn->size[1];
    emxEnsureCapacity_real_T(cn1, i);
    for (j = 0; j < n; j++) {
      ibcol = j * m;
      qEnd = j * iy;
      for (b_i = 0; b_i < m; b_i++) {
        cn1->data[ibcol + b_i] = 0.0;
      }

      for (k = 0; k < iy; k++) {
        vstride = k * m;
        t = cn->data[qEnd + k];
        for (b_i = 0; b_i < m; b_i++) {
          i = ibcol + b_i;
          cn1->data[i] += t * b_y->data[vstride + b_i];
        }
      }
    }
  }

  emxFree_real_T(&cn);

  // display(D);
  // display(V);
  t = 1.0 - pobs / 2.0;
  if ((t >= 0.0) && (t <= 1.0)) {
    t = -1.4142135623730951 * erfcinv(2.0 * t);
  } else {
    t = rtNaN;
  }

  i = b_y->size[0] * b_y->size[1];
  b_y->size[0] = cn1->size[0];
  b_y->size[1] = cn1->size[1];
  emxEnsureCapacity_real_T(b_y, i);
  qEnd = cn1->size[0] * cn1->size[1];
  for (i = 0; i < qEnd; i++) {
    b_y->data[i] = cn1->data[i] * cn1->data[i];
  }

  ibcol = b_y->size[1];
  if ((b_y->size[0] == 0) || (b_y->size[1] == 0)) {
    outsize_idx_0 = static_cast<unsigned int>(b_y->size[0]);
    i = a1->size[0];
    a1->size[0] = static_cast<int>(outsize_idx_0);
    emxEnsureCapacity_real_T(a1, i);
    qEnd = static_cast<int>(outsize_idx_0);
    for (i = 0; i < qEnd; i++) {
      a1->data[i] = 0.0;
    }
  } else {
    vstride = b_y->size[0];
    i = a1->size[0];
    a1->size[0] = b_y->size[0];
    emxEnsureCapacity_real_T(a1, i);
    for (j = 0; j < vstride; j++) {
      a1->data[j] = b_y->data[j];
    }

    for (k = 2; k <= ibcol; k++) {
      iy = (k - 1) * vstride;
      for (j = 0; j < vstride; j++) {
        a1->data[j] += b_y->data[iy + j];
      }
    }
  }

  emxFree_real_T(&b_y);
  iy = a1->size[0];
  for (k = 0; k < iy; k++) {
    a1->data[k] = std::sqrt(a1->data[k]);
  }

  qEnd = cn1->size[0];
  i = a1->size[0];
  a1->size[0] = cn1->size[0];
  emxEnsureCapacity_real_T(a1, i);
  for (i = 0; i < qEnd; i++) {
    a1->data[i] *= -t;
  }

  i = b1->size[0];
  b1->size[0] = a1->size[0];
  emxEnsureCapacity_real_T(b1, i);
  qEnd = a1->size[0];
  for (i = 0; i < qEnd; i++) {
    b1->data[i] = -a1->data[i];
  }

  // [p1,e1]=qscmvnv( 5000, r1, a1, cn1, b1);
  // [pvalue,evalue]=qscmvnvR(SimNumber, r1, a1, cn1, b1);
}

//
// File trailer for cnmatrix.cpp
//
// [EOF]
//
