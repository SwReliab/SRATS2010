// numlib.h

#ifndef NUMLIB_H
#define NUMLIB_H

#define INT_ZERO (1.0e-7)
#define INT_EPS (1.0e-6)
#define MTTF_MAXITE 20

#ifdef __cplusplus
extern "C" {
#endif

  /* gamma func */

  double loggamma(double x);
  double gamma(double x);
  double psi(double x);
  double polygamma(int n, double x);

  double fact(int s);
  double logfact(int s);

  /* gamma dist */

  double p_gamma(double a, double x, double loggamma_a);
  double q_gamma(double a, double x, double loggamma_a);

  /* Normal dist */

  double p_normal(double x);
  double q_normal(double x);
  double d_normal(double x);

  /* findquantile */

  double findQuantile(double q,
		      double (*cdf)(double, int, double*),
		      int npara, double *para);
  
#ifdef __cplusplus
}
#endif

#endif

