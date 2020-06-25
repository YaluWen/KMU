//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: mtimes.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 14-Apr-2020 21:40:05
//

// Include Files
#include "mtimes.h"
#include "qscmvnvR.h"
#include "qscmvnvR_emxutil.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : const emxArray_creal_T *A
//                const emxArray_creal_T *B
//                emxArray_creal_T *C
// Return Type  : void
//
void mtimes(const emxArray_creal_T *A, const emxArray_creal_T *B,
            emxArray_creal_T *C)
{
  int m;
  int inner;
  int i;
  int k;
  int aoffset;
  int B_re_tmp;
  m = A->size[0];
  inner = A->size[1];
  i = C->size[0];
  C->size[0] = A->size[0];
  emxEnsureCapacity_creal_T(C, i);
  for (i = 0; i < m; i++) {
    C->data[i].re = 0.0;
    C->data[i].im = 0.0;
  }

  for (k = 0; k < inner; k++) {
    aoffset = k * m;
    for (i = 0; i < m; i++) {
      B_re_tmp = aoffset + i;
      C->data[i].re += B->data[k].re * A->data[B_re_tmp].re - B->data[k].im *
        A->data[B_re_tmp].im;
      C->data[i].im += B->data[k].re * A->data[B_re_tmp].im + B->data[k].im *
        A->data[B_re_tmp].re;
    }
  }
}

//
// File trailer for mtimes.cpp
//
// [EOF]
//
