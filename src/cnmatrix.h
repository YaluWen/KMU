//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: cnmatrix.h
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:49:06
//
#ifndef CNMATRIX_H
#define CNMATRIX_H

// Include Files
#include <cstddef>
#include <cstdlib>
#include "rtwtypes.h"
//#include "cnmatrix_types.h"
#include "qscmvnvR_types.h"


// Function Declarations
extern void cnmatrix(double NoKernel, double pobs, const emxArray_real_T *V,
                     const emxArray_real_T *D, emxArray_real_T *r1,
                     emxArray_real_T *a1, emxArray_real_T *cn1, emxArray_real_T *
                     b1);

#endif

//
// File trailer for cnmatrix.h
//
// [EOF]
//
