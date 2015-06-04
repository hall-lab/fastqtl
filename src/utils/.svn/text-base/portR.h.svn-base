//
//  portR.h
//  portR
//
//  Created by Halit Ongen on 10/05/14.
//  Copyright (c) 2014 Halit Ongen. All rights reserved.
//

#ifndef portR_portR_h
#define portR_portR_h
#include <limits.h>
#include <math.h>
#include <float.h>

#define ML_POSINF       (1.0 / 0.0)
#define ML_NEGINF       ((-1.0) / 0.0)
#define ML_NAN          (0.0 / 0.0)
#ifndef M_LN2
#define M_LN2		0.693147180559945309417232121458	/* ln(2) */
#endif
#ifndef M_LOG10_2
#define M_LOG10_2	0.301029995663981195213738894724
#endif
#ifndef M_SQRT_PI
#define M_SQRT_PI	1.772453850905516027298167483341	/* sqrt(pi) */
#endif
#ifndef M_LN_SQRT_2PI
#define M_LN_SQRT_2PI	0.918938533204672741780329736406	/* log(sqrt(2*pi)) */
#endif
#define R_D__0	(log_p ? ML_NEGINF : 0.)		/* 0 */
#define R_D__1	(log_p ? 0. : 1.)			/* 1 */
#define R_D_exp(x)	(log_p	?  (x)	 : exp(x))
#define R_DT_0	(lower_tail ? R_D__0 : R_D__1)		/* 0 */
#define R_DT_1	(lower_tail ? R_D__1 : R_D__0)		/* 1 */
#define ISNAN(x) (isnan(x)!=0)
#define ML_VALID(x)	(!ISNAN(x))
#define R_P_bounds_01(x, x_min, x_max)  \
    if(x <= x_min) return R_DT_0;       \
    if(x >= x_max) return R_DT_1

#undef min
#define min(a,b) ((a < b)?a:b)
#undef max
#define max(a,b) ((a > b)?a:b)

#undef R_Log1_Exp
#define R_Log1_Exp(x)   ((x) > -M_LN2 ? log(-rexpm1(x)) : log1p(-exp(x)))

namespace portR {
    double pf(double x, double df1, double df2, int lower_tail, int log_p);
    double pbeta(double x, double pin, double qin, int lower_tail, int log_p);
    double pbeta_raw(double x, double pin, double qin, int lower_tail, int log_p);
    int Rf_i1mach(int);
    double Rf_d1mach(int);
    double logspace_add (double logx, double logy);
    double fmax2(double x, double y);
    
    double bfrac(double, double, double, double, double, double, int log_p);
    void bgrat(double, double, double, double, double *, double, int *, bool log_w);
    double grat_r(double a, double x, double r, double eps);
    double apser(double, double, double, double);
    double bpser(double, double, double, double, int log_p);
    double basym(double, double, double, double, int log_p);
    double fpser(double, double, double, double, int log_p);
    double bup(double, double, double, double, int, double, int give_log);
    double exparg(int);
    double psi(double);
    double gam1(double);
    double gamln1(double);
    double betaln(double, double);
    double algdiv(double, double);
    double brcmp1(int, double, double, double, double, int give_log);
    double brcomp(double, double, double, double, int log_p);
    double rlog1(double);
    double bcorr(double, double);
    double gamln(double);
    double alnrel(double);
    double esum(int, double, int give_log);
    double erf__(double);
    double rexpm1(double);
    double erfc1(int, double);
    double gsumln(double, double);
}
#endif
