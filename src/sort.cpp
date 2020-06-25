//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sort.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 14-Apr-2020 21:40:05
//

// Include Files
#include "sort.h"
#include "qscmvnvR.h"
#include "qscmvnvR_emxutil.h"
#include "rt_nonfinite.h"
#include "sortLE.h"

// Function Definitions

//
// Arguments    : emxArray_creal_T *x
//                emxArray_int32_T *idx
// Return Type  : void
//
void sort(emxArray_creal_T *x, emxArray_int32_T *idx)
{
  int n;
  int i2;
  int i;
  emxArray_int32_T *iwork;
  int k;
  emxArray_creal_T *xwork;
  int j;
  int pEnd;
  int p;
  int q;
  int qEnd;
  int kEnd;
  n = x->size[1] + 1;
  i2 = idx->size[0] * idx->size[1];
  idx->size[0] = 1;
  idx->size[1] = x->size[1];
  emxEnsureCapacity_int32_T(idx, i2);
  i = x->size[1];
  for (i2 = 0; i2 < i; i2++) {
    idx->data[i2] = 0;
  }

  if (x->size[1] != 0) {
    emxInit_int32_T(&iwork, 1);
    i2 = iwork->size[0];
    iwork->size[0] = x->size[1];
    emxEnsureCapacity_int32_T(iwork, i2);
    i2 = x->size[1] - 1;
    for (k = 1; k <= i2; k += 2) {
      if (sortLE(x, k, k + 1)) {
        idx->data[k - 1] = k;
        idx->data[k] = k + 1;
      } else {
        idx->data[k - 1] = k + 1;
        idx->data[k] = k;
      }
    }

    if ((x->size[1] & 1) != 0) {
      idx->data[x->size[1] - 1] = x->size[1];
    }

    i = 2;
    while (i < n - 1) {
      i2 = i << 1;
      j = 1;
      for (pEnd = i + 1; pEnd < n; pEnd = qEnd + i) {
        p = j;
        q = pEnd;
        qEnd = j + i2;
        if (qEnd > n) {
          qEnd = n;
        }

        k = 0;
        kEnd = qEnd - j;
        while (k + 1 <= kEnd) {
          if (sortLE(x, idx->data[p - 1], idx->data[q - 1])) {
            iwork->data[k] = idx->data[p - 1];
            p++;
            if (p == pEnd) {
              while (q < qEnd) {
                k++;
                iwork->data[k] = idx->data[q - 1];
                q++;
              }
            }
          } else {
            iwork->data[k] = idx->data[q - 1];
            q++;
            if (q == qEnd) {
              while (p < pEnd) {
                k++;
                iwork->data[k] = idx->data[p - 1];
                p++;
              }
            }
          }

          k++;
        }

        for (k = 0; k < kEnd; k++) {
          idx->data[(j + k) - 1] = iwork->data[k];
        }

        j = qEnd;
      }

      i = i2;
    }

    emxFree_int32_T(&iwork);
    emxInit_creal_T(&xwork, 1);
    i2 = xwork->size[0];
    xwork->size[0] = x->size[1];
    emxEnsureCapacity_creal_T(xwork, i2);
    for (k = 0; k <= n - 2; k++) {
      xwork->data[k] = x->data[k];
    }

    for (k = 0; k <= n - 2; k++) {
      x->data[k] = xwork->data[idx->data[k] - 1];
    }

    emxFree_creal_T(&xwork);
  }
}

//
// File trailer for sort.cpp
//
// [EOF]
//
