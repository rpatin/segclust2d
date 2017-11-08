#include <R.h>
#include <Rinternals.h>
#include <stdlib.h> // for NULL
#include <R_ext/Rdynload.h>

/* FIXME:
Check these declarations against the C/Fortran source code.
*/

/* .Call calls */
extern SEXP _segtools_apply_rowSums(SEXP, SEXP);
extern SEXP _segtools_arma_repmat(SEXP, SEXP, SEXP);
extern SEXP _segtools_colsums_sapply(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _segtools_DynProg_algo_cpp(SEXP, SEXP);
extern SEXP _segtools_Gmixt_algo_cpp(SEXP, SEXP, SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _segtools_Gmixt_simultanee_fullcpp(SEXP, SEXP, SEXP, SEXP, SEXP);
extern SEXP _segtools_logdens_simultanee_cpp(SEXP, SEXP, SEXP, SEXP);

static const R_CallMethodDef CallEntries[] = {
  {"_segtools_apply_rowSums",            (DL_FUNC) &_segtools_apply_rowSums,            2},
  {"_segtools_arma_repmat",              (DL_FUNC) &_segtools_arma_repmat,              3},
  {"_segtools_colsums_sapply",           (DL_FUNC) &_segtools_colsums_sapply,           5},
  {"_segtools_DynProg_algo_cpp",         (DL_FUNC) &_segtools_DynProg_algo_cpp,         2},
  {"_segtools_Gmixt_algo_cpp",           (DL_FUNC) &_segtools_Gmixt_algo_cpp,           7},
  {"_segtools_Gmixt_simultanee_fullcpp", (DL_FUNC) &_segtools_Gmixt_simultanee_fullcpp, 5},
  {"_segtools_logdens_simultanee_cpp",   (DL_FUNC) &_segtools_logdens_simultanee_cpp,   4},
  {NULL, NULL, 0}
};

void R_init_segtools(DllInfo *dll)
{
  R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
  R_useDynamicSymbols(dll, FALSE);
}
