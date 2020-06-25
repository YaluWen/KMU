//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: sum.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 14-Apr-2020 21:40:05
//

// Include Files
#include "sum.h"
#include "qscmvnvR.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : const emxArray_creal_T *x
// Return Type  : creal_T
//
creal_T sum(const emxArray_creal_T *x)
{
  creal_T y;
  int vlen;
  int k;
  vlen = x->size[1];
  if (x->size[1] == 0) {
    y.re = 0.0;
    y.im = 0.0;
  } else {
    y = x->data[0];
    for (k = 2; k <= vlen; k++) {
      y.re += x->data[k - 1].re;
      y.im += x->data[k - 1].im;
    }
  }

  return y;
}

//
// File trailer for sum.cpp
//
// [EOF]
//
