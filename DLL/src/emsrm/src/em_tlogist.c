#include <stdio.h>
#include <math.h>
#ifndef _MSC_VER
#include "numlib.h"
#else
#include "../include/numlib.h"
#endif

/////////////////////////////////////////////
// EM method
/////////////////////////////////////////////

double tlogist_pdf(double t, double shape, double scale) {
  double y = exp(-(t-scale)/shape);
  return (exp(-t/shape)+y)/(shape*(1.0+y)*(1.0+y));
}

double tlogist_cdf(double t, double shape, double scale) {
  return (1.0-exp(-t/shape))/(1.0+exp(-(t-scale)/shape));
}

double tlogist_mvf(double t, double omega, double shape, double scale) {
  return omega * tlogist_cdf(t, shape, scale);
}

double tlogist_cdf_v(double t, int npara, double *para) {
	return tlogist_cdf(t, para[1], para[2]);
}

void tlogist_inverse_mvf
(double *value, double *time, int *npara, double *para) {
  *time = findQuantile(*value/para[0], tlogist_cdf_v, *npara, para);
}

void tlogist_rate_series
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
    lambda[i] = omega * tlogist_pdf(t, shape, scale);
  }
}

double tlogist_creli(double x, double s, 
		   double omega, double shape, double scale,
		   double sval, double ffp) {
  double tmp;
  tmp = tlogist_mvf(x + s, omega, shape, scale) - sval;
  tmp = (exp(-tmp) - ffp) / (1.0 - ffp);
  return tmp;
}

void tlogist_mttf
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
  
  sval = tlogist_mvf(*stime, omega, shape, scale);
  ffp = exp(-(omega - sval));

  rev = sval - log(*r + (1.0-*r)*ffp);
  tlogist_inverse_mvf(&rev, &b, npara, para);
  b = b - *stime;
  
  // trapezoidal inte
  h = b - a;
  t = h * (tlogist_creli(a, *stime, omega, shape, scale, sval, ffp)
	   + tlogist_creli(b, *stime, omega, shape, scale, sval, ffp)) / 2.0;

  for (k = 1; k <= MTTF_MAXITE; k++) {
    n = 2 * n;
    h = h / 2.0;
    s = 0.0;
    for (i=1; i<=n-1; i+=2) {
      s += tlogist_creli(a + i * h, *stime, omega, shape, scale, sval, ffp);
    }
    tn = t/2.0 + h * s;
    if (fabs(tn - t) < INT_EPS*fabs(t)) {
      break;
    }
    t = tn;
  }
  *ss = tn;
}

void tlogist_mvf_series
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
    mean[i] = tlogist_mvf(t, omega, shape, scale);
  }
}

void em_tlogist_emstep
(int *dsize, double *time, double *num, int *type, 
 int *npara, double *para, double *pdiff, double *retllf, double *total) {

  int j;

  double x, t, y, y0;
  double tmp1, tmp2, tmp3;
  double llf;
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
  y0 = exp(-scale/shape);
  omega = omega * (1.0 + y0);

  en1 = 0.0;
  en2 = 0.0;
  en3 = 0.0;
  llf = 0.0;
  g00 = 1.0/(1.0+y0);
  g01 = 1.0/(2.0*(1.0+y0)*(1.0+y0));
  g02 = (1.0+(1.0+log(y0))*y0)/((1.0+y0)*(1.0+y0));
  t = time[0];
  x = num[0];
  y = exp((t-scale)/shape);
  g10 = g00;
  g11 = g01;
  g12 = g02;
  g20 = 1.0/(1.0+y);
  g21 = 1.0/(2.0*(1.0+y)*(1.0+y));
  g22 = (1.0+(1.0+log(y))*y)/((1.0+y)*(1.0+y));
  if (x != 0.0) {
    tmp1 = g10 - g20;
    tmp2 = g11 - g21;
    tmp3 = g12 - g22;
    en1 += x;
    en2 += x * tmp2 / tmp1;
    en3 += x * tmp3 / tmp1;
    llf += x * (log(tmp1) - log(g00)) - loggamma(x+1);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += 1.0/(1.0+y);
    en3 += (y-1.0)*log(y)/(1.0+y);
    llf += log(tlogist_pdf(t, shape, scale));
  }
  for (j=1; j<*dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      y = exp((t-scale)/shape);
      g10 = g20;
      g11 = g21;
      g12 = g22;
      g20 = 1.0/(1.0+y);
      g21 = 1.0/(2.0*(1.0+y)*(1.0+y));
      g22 = (1.0+(1.0+log(y))*y)/((1.0+y)*(1.0+y));
    }
    if (x != 0.0) {
      tmp1 = g10 - g20;
      tmp2 = g11 - g21;
      tmp3 = g12 - g22;
      en1 += x;
      en2 += x * tmp2 / tmp1;
      en3 += x * tmp3 / tmp1;
      llf += x * (log(tmp1) - log(g00))- loggamma(x+1);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += 1.0/(1.0+y);
      en3 += (y-1.0)*log(y)/(1.0+y);
      llf += log(tlogist_pdf(t, shape, scale));
    }
  }
  llf += (log(omega) + log(g00)) * en1;  // en1 is total number of faults
  *total = en1 + omega * g20;
  en1 += omega * (1.0 - g00 + g20);  // g00 is the first, g10 is the last
  en2 += omega * (0.5 - g01 + g21);  // g01 is the first, g11 is the last
  en3 += omega * (1.0 - g02 + g22);  // g02 is the first, g12 is the last
  llf += - omega * (g00 - g20);
		
  // M-step
  shape = shape*en3/en1;
  scale = scale + shape*(log(en1/2.0) - log(en2));
  y0 = exp(-scale/shape);
  omega = en1 / (1.0 + y0);
  // *total = omega;

  pdiff[0] = omega - para[0];
  pdiff[1] = shape - para[1];
  pdiff[2] = scale - para[2];

  para[0] = omega;
  para[1] = shape;
  para[2] = scale;

  *retllf = llf;
}

