//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: relop.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 14-Apr-2020 21:40:05
//

// Include Files
#include "relop.h"
#include "qscmvnvR.h"
#include "qscmvnvR_rtwutil.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <math.h>

// Function Definitions

//
// Arguments    : double x
//                double y
// Return Type  : boolean_T
//
boolean_T iseq(double x, double y)
{
  boolean_T p;
  double absx;
  int exponent;
  absx = std::abs(y / 2.0);
  if ((!rtIsInf(absx)) && (!rtIsNaN(absx))) {
    if (absx <= 2.2250738585072014E-308) {
      absx = 4.94065645841247E-324;
    } else {
      frexp(absx, &exponent);
      absx = std::ldexp(1.0, exponent - 53);
    }
  } else {
    absx = rtNaN;
  }

  if ((std::abs(y - x) < absx) || (rtIsInf(x) && rtIsInf(y) && ((x > 0.0) == (y >
         0.0)))) {
    p = true;
  } else {
    p = false;
  }

  return p;
}

//
// Arguments    : const creal_T a
//                double b
// Return Type  : boolean_T
//
boolean_T relop(const creal_T a, double b)
{
  boolean_T p;
  double ma;
  double absai;
  double absbr;
  double Ma;
  if (rtIsNaN(b)) {
    p = false;
  } else if (rtIsNaN(a.re) || rtIsNaN(a.im)) {
    p = true;
  } else {
    ma = std::abs(a.re);
    if ((ma > 8.9884656743115785E+307) || (std::abs(a.im) >
         8.9884656743115785E+307)) {
      absai = rt_hypotd_snf(a.re / 2.0, a.im / 2.0);
      absbr = std::abs(b) / 2.0;
    } else {
      absai = rt_hypotd_snf(a.re, a.im);
      absbr = std::abs(b);
    }

    if (absai == absbr) {
      absai = std::abs(a.im);
      absbr = std::abs(b);
      if (ma > absai) {
        Ma = ma;
        ma = absai;
      } else {
        Ma = absai;
      }

      if (absbr > 0.0) {
        absai = absbr;
        absbr = 0.0;
      } else {
        absai = 0.0;
      }

      if (Ma > absai) {
        ma = Ma;
        absbr = absai;
      } else {
        if (Ma < absai) {
          if (ma > absbr) {
            absbr = absai - Ma;
            ma *= ma / 2.0 / (Ma / 2.0 + absai / 2.0);
          } else {
            ma = Ma;
            absbr = absai;
          }
        }
      }

      absai = ma;
      if (ma == absbr) {
        absai = rt_atan2d_snf(a.im, a.re);
        absbr = rt_atan2d_snf(0.0, b);
        if (absai == absbr) {
          if (absai > 0.78539816339744828) {
            if (absai > 2.3561944901923448) {
              absai = -a.im;
              absbr = -0.0;
            } else {
              absai = -a.re;
              absbr = -b;
            }
          } else if (absai > -0.78539816339744828) {
            absai = a.im;
            absbr = 0.0;
          } else if (absai > -2.3561944901923448) {
            absai = a.re;
            absbr = b;
          } else {
            absai = -a.im;
            absbr = -0.0;
          }

          if (absai == absbr) {
            absai = 0.0;
            absbr = 0.0;
          }
        }
      }
    }

    p = (absai < absbr);
  }

  return p;
}

//
// File trailer for relop.cpp
//
// [EOF]
//
