//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: eml_primes_core.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 14-Apr-2020 21:40:05
//

// Include Files
#include "eml_primes_core.h"
#include "qscmvnvR.h"
#include "rt_nonfinite.h"

// Function Definitions

//
// Arguments    : unsigned int a
// Return Type  : unsigned int
//
unsigned int intsqrt(unsigned int a)
{
  unsigned int r;
  unsigned int rhi;
  unsigned int rlo;
  boolean_T exitg1;
  unsigned int rr;
  rhi = (a >> 1) + 1U;
  if (65535U < rhi) {
    rhi = 65535U;
  }

  if (a >= rhi * rhi) {
    r = rhi;
  } else {
    rlo = 0U;
    r = rhi >> 1;
    exitg1 = false;
    while ((!exitg1) && (r > rlo)) {
      rr = r * r;
      if (rr == a) {
        exitg1 = true;
      } else {
        if (rr > a) {
          rhi = r;
        } else {
          rlo = r;
        }

        r = rlo + ((rhi - rlo) >> 1);
        if (r < rlo) {
          r = MAX_uint32_T;
        }
      }
    }
  }

  return r;
}

//
// File trailer for eml_primes_core.cpp
//
// [EOF]
//
