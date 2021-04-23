#include <R.h>
#include <Rinternals.h>

#include "robustarima.h"

#include <R_ext/Rdynload.h>

/* Fortran interface descriptions: */

static R_NativePrimitiveArgType s_regafe_t[49] = {
  REALSXP,  /* x, */
  REALSXP,  /* y, */
  INTSXP,   /* n, */
  INTSXP,   /* m, */
  INTSXP ,  /* idif, */
  INTSXP ,  /* isp, */
  INTSXP ,  /* nsd, */
  INTSXP ,  /* ip, */
  INTSXP ,  /* indar, */
  INTSXP ,  /* ipfin, */
  INTSXP ,  /* iqfin, */
  INTSXP ,  /* indth, */
  INTSXP ,  /* interc, */
  INTSXP ,  /* kopt, */
  REALSXP,  /* phiopt, */
  REALSXP,  /* theta, */
  REALSXP,  /* thetas, */
  REALSXP,  /* betaopt, */
  REALSXP,  /* sigmau, */
  REALSXP,  /* bcov, */
  REALSXP,  /* zcor, */
  REALSXP,  /* zcor1, */
  REALSXP,  /* sigmadif, */
  REALSXP,  /* cck, */
  REALSXP,  /* sigfil, */
  REALSXP,  /* xy, */
  REALSXP,  /* yhat, */
  REALSXP,  /* uhat, */
  REALSXP,  /* epshat, */
  REALSXP,  /* st, */
  REALSXP,  /* epspred, */
  INTSXP,   /* npred, */
  REALSXP,  /* tauef, */
  INTSXP,   /* infnew, */
  INTSXP,   /* ndim1, */
  INTSXP,   /* ndim2, */
  REALSXP,  /* work1, */
  INTSXP,   /* nw1, */
  INTSXP,   /* iwork1, */
  INTSXP,   /* niw1, */
  REALSXP,  /* work2, */
  INTSXP,   /* nw2, */
  INTSXP,   /* iwork2, */
  INTSXP,   /* niw2, */
  REALSXP,  /* utol, */
  INTSXP,   /* maxfev, */
  REALSXP,  /* epsmch, */
  REALSXP,  /* dwarf, */
  INTSXP    /* n0 */
} ;

static R_NativePrimitiveArgType s_outlfe_t[33] = {
  REALSXP,  /* x, */
  REALSXP,  /* y, */
  INTSXP,   /* n, */
  INTSXP,   /* m, */
  INTSXP,   /* idif, */
  INTSXP,   /* isp, */
  INTSXP,   /* nsd, */
  INTSXP,   /* k, */
  INTSXP,   /* q, */
  INTSXP,   /* indth, */
  REALSXP,  /* beta, */
  REALSXP,  /* phidif, */
  REALSXP,  /* theta, */
  REALSXP,  /* thetas, */
  REALSXP,  /* sigmadif, */
  INTSXP,   /* indio, */
  REALSXP,  /* cck, */
  REALSXP,  /* sigfil, */
  REALSXP,  /* critv, */
  INTSXP,   /* nout, */
  INTSXP,   /* indtipo, */
  INTSXP,   /* t0, */
  REALSXP,  /* wout, */
  REALSXP,  /* lambda, */
  REALSXP,  /* sigmau0, */
  REALSXP,  /* sigmau, */
  INTSXP,   /* idim, */
  REALSXP,  /* work, */
  INTSXP,   /* idimw, */
  INTSXP,   /* iwork, */
  INTSXP,   /* idimiw, */
  INTSXP,   /* ierror, */
  INTSXP    /* n0 */
} ;

static R_FortranMethodDef fortranMethods[] = {
  {"s_regafe", (DL_FUNC) &F77_SUB(s_regafe), 49, s_regafe_t},
  {"s_outlfe", (DL_FUNC) &F77_SUB(s_outlfe), 33, s_outlfe_t},
  {NULL, NULL, 0}
} ;

void
R_init_robustarima(DllInfo *dll)
{
  R_registerRoutines(dll,
    NULL /*cMethods*/, NULL /*callMethods*/,
    fortranMethods, NULL /*externalMethods*/) ;
  R_useDynamicSymbols(dll, FALSE) ;
}
