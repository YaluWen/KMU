//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: cnmatrix_initialize.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:49:06
//

// Include Files
#include "cnmatrix_initialize.h"
#include "cnmatrix.h"
#include "cnmatrix_data.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void cnmatrix_initialize()
{
  rt_InitInfAndNaN();
  isInitialized_cnmatrix = true;
}

//
// File trailer for cnmatrix_initialize.cpp
//
// [EOF]
//
