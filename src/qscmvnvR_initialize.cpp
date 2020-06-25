//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: qscmvnvR_initialize.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 14-Apr-2020 21:40:05
//

// Include Files
#include "qscmvnvR_initialize.h"
#include "eml_rand_mt19937ar_stateful.h"
#include "qscmvnvR.h"
#include "qscmvnvR_data.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : void
// Return Type  : void
//
void qscmvnvR_initialize()
{
  rt_InitInfAndNaN();
  c_eml_rand_mt19937ar_stateful_i();
  isInitialized_qscmvnvR = true;
}

//
// File trailer for qscmvnvR_initialize.cpp
//
// [EOF]
//
