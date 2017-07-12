#include <math.h>
#ifndef _MSC_VER
#include "numlib.h"
#else
#include "../include/numlib.h"
#endif

/////////////////////////////////////////////
// EM method
/////////////////////////////////////////////

double tnorm_pdf(double t, double shape, double scale) {
  return d_normal((t-scale)/shape) / q_normal(-scale/shape)/shape;
}

double tnorm_cdf(double t, double shape, double scale) {
  return 1.0 - q_normal((t-scale)/shape)/q_normal(-scale/shape);
}

double tnorm_mvf(double t, double omega, double shape, double scale) {
  return omega * tnorm_cdf(t, shape, scale);
}

double tnorm_cdf_v(double t, int npara, double *para) {
	return tnorm_cdf(t, para[1], para[2]);
}

void tnorm_inverse_mvf
(double *value, double *time, int *npara, double *para) {
  *time = findQuantile(*value/para[0], tnorm_cdf_v, *npara, para);
}

void tnorm_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para) {
  int i;

  double t;
  double omega, shape, scale;

  if (*npara != 3) {
    return;
  }
  omega = para[0];
  shape = para[1];
  scale = para[2];

  t = 0;
  for (i=0; i<*dsize; i++) {
    t += time[i];
    lambda[i] = omega * tnorm_pdf(t, shape, scale);
  }
}

double tnorm_creli(double x, double s, 
		   double omega, double shape, double scale,
		   double sval, double ffp) {
  double tmp;
  tmp = tnorm_mvf(x + s, omega, shape, scale) - sval;
  tmp = (exp(-tmp) - ffp) / (1.0 - ffp);
  return tmp;
}

void tnorm_mttf
(double *stime, double *r, double *ss, int *npara, double *para) {
  int i, k;
  int n = 1;
  double s, h, t, tn;
  double ffp, sval;
  double rev, b, a = 0.0;
  double omega, shape, scale;

  if (*npara != 3) {
    *ss = 0.0;
    return;
  }
  omega = para[0];
  shape = para[1];
  scale = para[2];
  
  sval = tnorm_mvf(*stime, omega, shape, scale);
  ffp = exp(-(omega - sval));

  rev = sval - log(*r + (1.0-*r)*ffp);
  tnorm_inverse_mvf(&rev, &b, npara, para);
  b = b - *stime;
  
  // trapezoidal inte
  h = b - a;
  t = h * (tnorm_creli(a, *stime, omega, shape, scale, sval, ffp)
	   + tnorm_creli(b, *stime, omega, shape, scale, sval, ffp)) / 2.0;

  for (k = 1; k <= MTTF_MAXITE; k++) {
    n = 2 * n;
    h = h / 2.0;
    s = 0.0;
    for (i=1; i<=n-1; i+=2) {
      s += tnorm_creli(a + i * h, *stime, omega, shape, scale, sval, ffp);
    }
    tn = t/2.0 + h * s;
    if (fabs(tn - t) < INT_EPS*fabs(t)) {
      break;
    }
    t = tn;
  }
  *ss = tn;
}

void tnorm_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para) {
  int i;

  double t;
  double omega;
  double shape;
  double scale;

  if (*npara != 3) {
    return;
  }
  omega = para[0];
  shape = para[1];
  scale = para[2];

  t = 0;
  for (i=0; i<*dsize; i++) {
    t += time[i];
    mean[i] = tnorm_mvf(t, omega, shape, scale);
  }
}

void em_tnorm_emstep
(int *dsize, double *time, double *num, int *type, 
 int *npara, double *para, double *pdiff, double *retllf, double *total) {

  int j;
  double x, t, y;
  double tmp1, tmp2, tmp3;
  double tmp, llf;
  double g00, g01, g02;
  double g10, g11, g12;
  double g20, g21, g22;
  double en1, en2, en3;
		
  double omega;
  double shape;
  double scale;

  if (*npara != 3) {
    *retllf = 0.0;
    return;
  }
  omega = para[0];
  shape = para[1];
  scale = para[2];

  // E-step
  omega = omega / q_normal(-scale/shape);

  en1 = 0.0;
  en2 = 0.0;
  en3 = 0.0;
  llf = 0.0;
  tmp = d_normal(-scale/shape);
  g00 = q_normal(-scale/shape);
  g01 = shape*tmp + scale*g00;
  g02 = (scale*shape)*tmp + (shape*shape + scale*scale)*g00;
  t = time[0];
  x = num[0];
  y = (t-scale)/shape;
  tmp = d_normal(y);
  g10 = g00;
  g11 = g01;
  g12 = g02;
  g20 = q_normal(y);
  g21 = shape*tmp + scale*g20;
  g22 = (shape*t + scale*shape)*tmp + (shape*shape + scale*scale)*g20;
  if (x != 0.0) {
    tmp1 = g10 - g20;
    tmp2 = g11 - g21;
    tmp3 = g12 - g22;
    en1 += x;
    en2 += x * tmp2 / tmp1;
    en3 += x * tmp3 / tmp1;
    llf += x * (log(tmp1) - log(g00)) - loggamma(x+1.0);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += t;
    en3 += t * t;
    llf += log(tnorm_pdf(t, shape, scale));
  }
  for (j=1; j<*dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      y = (t-scale)/shape;
      tmp = d_normal(y);
      g10 = g20;
      g11 = g21;
      g12 = g22;
      g20 = q_normal(y);
      g21 = shape*tmp + scale*g20;
      g22 = (shape*t + scale*shape)*tmp + (shape*shape + scale*scale)*g20;
    }
    if (x != 0.0) {
      tmp1 = g10 - g20;
      tmp2 = g11 - g21;
      tmp3 = g12 - g22;
      en1 += x;
      en2 += x * tmp2 / tmp1;
      en3 += x * tmp3 / tmp1;
      llf += x * (log(tmp1) - log(g00))- loggamma(x+1.0);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += t;
      en3 += t * t;
      llf += log(tnorm_pdf(t, shape, scale));
    }
  }
  llf += (log(omega) + log(g00)) * en1;  // en1 is total number of faults
  *total = en1 + omega * g20;
  en1 += omega * (1.0 - g00 + g20);  // g00 is the first, g10 is the last
  en2 += omega * (scale - g01 + g21);  // g01 is the first, g11 is the last
  en3 += omega * (shape*shape + scale*scale 
		  - g02 + g22);  // g02 is the first, g22 is the last
  llf += - omega * (g00 - g20);
  
  // M-step
  scale = en2 / en1;
  shape = sqrt(en3 / en1 - scale*scale);
  omega =  en1 * q_normal(-scale/shape);
  //*total = omega;

  pdiff[0] = omega - para[0];
  pdiff[1] = shape - para[1];
  pdiff[2] = scale - para[2];

  para[0] = omega;
  para[1] = shape;
  para[2] = scale;

  *retllf = llf;
}

