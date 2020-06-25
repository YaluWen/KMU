//
// Academic License - for use in teaching, academic research, and meeting
// course requirements at degree granting institutions only.  Not for
// government, commercial, or other organizational use.
// File: erfcinv.cpp
//
// MATLAB Coder version            : 4.3
// C/C++ source code generated on  : 15-Apr-2020 13:49:06
//

// Include Files
#include "erfcinv.h"
#include "cnmatrix.h"
#include "cnmatrix_rtwutil.h"
#include "rt_nonfinite.h"
#include <cmath>
#include <math.h>

// Function Declarations
static double eml_erfcore(double x);

// Function Definitions

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
// Arguments    : double y
// Return Type  : double
//
double erfcinv(double y)
{
  double x;
  int nIterations;
  boolean_T guard1 = false;
  double z;
  int k;
  double absx;
  double P;
  double R;
  double S;
  int eint;
  nIterations = 1;
  if (rtIsNaN(y)) {
    x = rtNaN;
  } else if (y == 0.0) {
    x = rtInf;
  } else if (y == 2.0) {
    x = rtMinusInf;
  } else {
    guard1 = false;
    if (y > 1.7) {
      z = std::sqrt(-std::log((2.0 - y) / 2.0));
      x = -(((1.641345311 * z + 3.429567803) * z + -1.624906493) * z +
            -1.970840454) / ((1.6370678 * z + 3.5438892) * z + 1.0);
      guard1 = true;
    } else if (y < 0.3) {
      z = std::sqrt(0.69314718055994529 - std::log(y));
      if (y < 1.0947644252537633E-47) {
        if (y < 7.7532508072625747E-267) {
          nIterations = 3;
        } else {
          nIterations = 2;
        }
      }

      x = (((1.641345311 * z + 3.429567803) * z + -1.624906493) * z +
           -1.970840454) / ((1.6370678 * z + 3.5438892) * z + 1.0);
      for (k = 0; k <= nIterations; k++) {
        z = -(eml_erfcore(x) - y) / (1.1283791670955126 * std::exp(-x * x));
        x -= z / (x * z + 1.0);
      }
    } else {
      z = (1.0 - y) * (1.0 - y);
      x = (1.0 - y) * (((-0.140543331 * z + 0.914624893) * z + -1.645349621) * z
                       + 0.886226899) / ((((0.012229801 * z + -0.329097515) * z
        + 1.442710462) * z + -2.118377725) * z + 1.0);
      guard1 = true;
    }

    if (guard1) {
      if (y > 1.7) {
        for (k = 0; k < 2; k++) {
          z = (eml_erfcore(-x) - (2.0 - y)) / (1.1283791670955126 * std::exp(-x *
            x));
          x -= z / (x * z + 1.0);
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
          absx = std::abs(x);
          if (rtIsNaN(x)) {
            z = x;
          } else if (rtIsInf(x)) {
            if (x < 0.0) {
              z = -1.0;
            } else {
              z = 1.0;
            }
          } else if (absx < 0.84375) {
            if (absx < 3.7252902984619141E-9) {
              if (absx < 2.8480945388892178E-306) {
                z = 0.125 * (8.0 * x + 1.0270333367641007 * x);
              } else {
                z = x + 0.12837916709551259 * x;
              }
            } else {
              z = x * x;
              z = x + x * ((z * (z * (z * (z * -2.3763016656650163E-5 +
                -0.0057702702964894416) + -0.02848174957559851) +
                                 -0.3250421072470015) + 0.12837916709551256) /
                           (z * (z * (z * (z * (z * -3.9602282787753681E-6 +
                0.00013249473800432164) + 0.0050813062818757656) +
                0.0650222499887673) + 0.39791722395915535) + 1.0));
            }
          } else if (absx < 1.25) {
            P = (absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) *
              ((absx - 1.0) * ((absx - 1.0) * -0.0021663755948687908 +
                               0.035478304325618236) + -0.11089469428239668) +
              0.31834661990116175) + -0.37220787603570132) + 0.41485611868374833)
              + -0.0023621185607526594;
            z = (absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) * ((absx - 1.0) *
              ((absx - 1.0) * ((absx - 1.0) * 0.011984499846799107 +
                               0.013637083912029051) + 0.12617121980876164) +
              0.071828654414196266) + 0.540397917702171) + 0.10642088040084423)
              + 1.0;
            if (x >= 0.0) {
              z = P / z + 0.84506291151046753;
            } else {
              z = -0.84506291151046753 - P / z;
            }
          } else if (absx > 6.0) {
            if (x < 0.0) {
              z = -1.0;
            } else {
              z = 1.0;
            }
          } else {
            z = 1.0 / (absx * absx);
            if (absx < 2.8571434020996094) {
              R = z * (z * (z * (z * (z * (z * (z * -9.8143293441691455 +
                -81.2874355063066) + -184.60509290671104) + -162.39666946257347)
                                 + -62.375332450326006) + -10.558626225323291) +
                       -0.69385857270718176) + -0.0098649440348471482;
              S = z * (z * (z * (z * (z * (z * (z * (z * -0.0604244152148581 +
                6.5702497703192817) + 108.63500554177944) + 429.00814002756783)
                          + 645.38727173326788) + 434.56587747522923) +
                            137.65775414351904) + 19.651271667439257) + 1.0;
            } else {
              R = z * (z * (z * (z * (z * (z * -483.5191916086514 +
                -1025.0951316110772) + -637.56644336838963) +
                                 -160.63638485582192) + -17.757954917754752) +
                       -0.799283237680523) + -0.0098649429247001;
              S = z * (z * (z * (z * (z * (z * (z * -22.440952446585818 +
                474.52854120695537) + 2553.0504064331644) + 3199.8582195085955)
                                 + 1536.729586084437) + 325.79251299657392) +
                       30.338060743482458) + 1.0;
            }

            if (!rtIsNaN(absx)) {
              P = frexp(absx, &eint);
              nIterations = eint;
            } else {
              P = absx;
              nIterations = 0;
            }

            z = std::floor(P * 2.097152E+6) / 2.097152E+6 * rt_powd_snf(2.0,
              static_cast<double>(nIterations));
            z = std::exp(-z * z - 0.5625) * std::exp((z - absx) * (z + absx) + R
              / S) / absx;
            if (x < 0.0) {
              z--;
            } else {
              z = 1.0 - z;
            }
          }

          z = (z - (1.0 - y)) / (1.1283791670955126 * std::exp(-x * x));
          x -= z / (x * z + 1.0);
        }
      }
    }
  }

  return x;
}

//
// File trailer for erfcinv.cpp
//
// [EOF]
//
