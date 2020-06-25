// Include Files
//#include "main.h"
#include "qscmvnvR.h"
#include "qscmvnvR_emxAPI.h"
#include "qscmvnvR_terminate.h"
#include "rt_nonfinite.h"
#include "cnmatrix.h"
#include "cnmatrix_emxAPI.h"
#include "cnmatrix_terminate.h"
#include "RcppArmadillo.h"

using namespace Rcpp;

// Function Declarations
static emxArray_real_T *argInit_Unboundedx1_real_T();
static double argInit_real_T();
//static emxArray_real_T *c_argInit_UnboundedxUnbounded_r();
static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(arma::mat CovVarAll);

// Function Definitions
arma::mat eigenvalues(arma::mat A) {
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A);
  //eigvec.print("Eigen Vectors");
  //eigval.print("Eigen Vectors");
  arma::mat eigvalmat(eigvec.n_rows,eigvec.n_cols);
  eigvalmat.fill(0);
  eigvalmat.diag()=eigval;
  return eigvalmat;
}
arma::mat eigenvec(arma::mat A) {
  arma::vec eigval;
  arma::mat eigvec;
  arma::eig_sym(eigval, eigvec, A);
  //eigvec.print("Eigen Vectors");
  //eigval.print("Eigen Vectors");
  return eigvec;
}

//
// Arguments    : void
// Return Type  : emxArray_real_T *
//
static emxArray_real_T *c_argInit_UnboundedxUnbounded_r(arma::mat CovVarAll)
{
  emxArray_real_T *result;
  int idx0;
  int idx1;
  
  // Set the size of the array.
  // Change this size to the value that the application requires.
  result = emxCreate_real_T(CovVarAll.n_rows, CovVarAll.n_cols);
  
  // Loop over the array to initialize each element.
  for (idx0 = 0; idx0 < result->size[0U]; idx0++) {
    for (idx1 = 0; idx1 < result->size[1U]; idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      result->data[idx0 + result->size[0] * idx1] = CovVarAll(idx0,idx1);
    }
  }
  
  return result;
}

// [[Rcpp::export]]
double main_UStatsPvalueMain(double NoKernel, double pobs, arma::mat CovVarAll,double SimNumber)
{
  arma::mat Dmatrix;
  arma::mat Vmatrix;
  Vmatrix=eigenvec(CovVarAll);
  Dmatrix=eigenvalues(CovVarAll);

  
  emxArray_real_T *r1;
  emxArray_real_T *a1;
  emxArray_real_T *cn1;
  emxArray_real_T *b1;

  emxArray_real_T *V;
  emxArray_real_T *D;
  
  emxInitArray_real_T(&r1, 2);
  emxInitArray_real_T(&a1, 1);
  emxInitArray_real_T(&cn1, 2);
  emxInitArray_real_T(&b1, 1);

  double p;
  double e;
  
  // Initialize function input argument 'V'.
  V = c_argInit_UnboundedxUnbounded_r(Vmatrix);

  // Initialize function input argument 'D'.
  D = c_argInit_UnboundedxUnbounded_r(Dmatrix);

  // Call the entry-point 'cnmatrix'.
  cnmatrix(NoKernel, pobs, V, D, r1, a1, cn1, b1);
 /* for (int idx0 = 0; idx0 < cn1->size[0U]; idx0++) {
    for (int idx1 = 0; idx1 < cn1->size[1U]; idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      Rcout<< cn1->data[idx0 + cn1->size[0] * idx1] << ",";
    }
    Rcout<<"\n";
  }
  for (int idx0 = 0; idx0 < a1->size[0U]; idx0++) {
      Rcout<< a1->data[idx0];
    Rcout<<"\n";
  }
  Rcout<<"bl\n";
  
  for (int idx0 = 0; idx0 < b1->size[0U]; idx0++) {
    Rcout<< b1->data[idx0];
    Rcout<<"\n";
  }
  
  for (int idx0 = 0; idx0 < r1->size[0U]; idx0++) {
    for (int idx1 = 0; idx1 < r1->size[1U]; idx1++) {
      // Set the value of the array element.
      // Change this value to the value that the application requires.
      Rcout<< r1->data[idx0 + r1->size[0] * idx1] << "t";
    }
    Rcout<<"\n";
  }
  Rcout<<"Nokernel"<<NoKernel;*/
  qscmvnvR(SimNumber, r1, a1, cn1, b1, &p, &e);
  
  emxDestroyArray_real_T(b1);
  emxDestroyArray_real_T(cn1);
  emxDestroyArray_real_T(a1);
  emxDestroyArray_real_T(r1);
  emxDestroyArray_real_T(V);
  emxDestroyArray_real_T(D);
  
  return p;
}



// [EOF]
//
