#include <math.h>
#ifndef _MSC_VER
#include "numlib.h"
#else
#include "../include/numlib.h"
#endif

/////////////////////////////////////////////
// EM method
/////////////////////////////////////////////

double lnorm_pdf(double t, double shape, double scale) {
  return d_normal((log(t)-scale)/shape)/(t*shape);
}

double lnorm_cdf(double t, double shape, double scale) {
  return p_normal((log(t)-scale)/shape);
}

double lnorm_mvf(double t, double omega, double shape, double scale) {
  return omega * lnorm_cdf(t, shape, scale);
}

double lnorm_cdf_v(double t, int npara, double *para) {
	return lnorm_cdf(t, para[1], para[2]);
}

void lnorm_inverse_mvf
(double *value, double *time, int *npara, double *para) {
  *time = findQuantile(*value/para[0], lnorm_cdf_v, *npara, para);
}

void lnorm_rate_series
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
    lambda[i] = omega * lnorm_pdf(t, shape, scale);
  }
}

double lnorm_creli(double x, double s, 
		   double omega, double shape, double scale,
		   double sval, double ffp) {
  double tmp;
  tmp = lnorm_mvf(x + s, omega, shape, scale) - sval;
  tmp = (exp(-tmp) - ffp) / (1.0 - ffp);
  return tmp;
}

void lnorm_mttf
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
  
  sval = lnorm_mvf(*stime, omega, shape, scale);
  ffp = exp(-(omega - sval));

  rev = sval - log(*r + (1.0-*r)*ffp);
  lnorm_inverse_mvf(&rev, &b, npara, para);
  b = b - *stime;
  
  // trapezoidal inte
  h = b - a;
  t = h * (lnorm_creli(a, *stime, omega, shape, scale, sval, ffp)
	   + lnorm_creli(b, *stime, omega, shape, scale, sval, ffp)) / 2.0;

  for (k = 1; k <= MTTF_MAXITE; k++) {
    n = 2 * n;
    h = h / 2.0;
    s = 0.0;
    for (i=1; i<=n-1; i+=2) {
      s += lnorm_creli(a + i * h, *stime, omega, shape, scale, sval, ffp);
    }
    tn = t/2.0 + h * s;
    if (fabs(tn - t) < INT_EPS*fabs(t)) {
      break;
    }
    t = tn;
  }
  *ss = tn;
}

void lnorm_mvf_series
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
    mean[i] = lnorm_mvf(t, omega, shape, scale);
  }
}

void em_lnorm_emstep
(int *dsize, double *time, double *num, int *type, 
 int *npara, double *para, double *pdiff, double *retllf, double *total) {

  int j;
  double x, t, y;
  double tmp1, tmp2, tmp3;
  double tmp, llf;
  double g00, g01, g02;
  double g10, g11, g12;
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
  t = time[0];
  x = num[0];
  y = (log(t)-scale)/shape;
  en1 = 0.0;
  en2 = 0.0;
  en3 = 0.0;
  llf = 0.0;
  tmp = d_normal(y);
  g00 = 1.0;
  g01 = scale;
  g02 = shape*shape + scale*scale;
  g10 = q_normal(y);
  g11 = shape*tmp + scale*g10;
  g12 = (shape*log(t) + scale*shape)*tmp 
    + (shape*shape + scale*scale)*g10;
  if (x != 0.0) {
    tmp1 = g00 - g10;
    tmp2 = g01 - g11;
    tmp3 = g02 - g12;
    en1 += x;
    en2 += x * tmp2 / tmp1;
    en3 += x * tmp3 / tmp1;
    llf += x * log(tmp1) - loggamma(x+1);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += log(t);
    en3 += log(t) * log(t);
    llf += log(lnorm_pdf(t, shape, scale));
  }
  for (j=1; j<*dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      y = (log(t)-scale)/shape;
      tmp = d_normal(y);
      g00 = g10;
      g01 = g11;
      g02 = g12;
      g10 = q_normal(y);
      g11 = shape*tmp + scale*g10;
      g12 = (shape*log(t) + scale*shape)*tmp 
	+ (shape*shape + scale*scale)*g10;
    }
    if (x != 0.0) {
      tmp1 = g00 - g10;
      tmp2 = g01 - g11;
      tmp3 = g02 - g12;
      en1 += x;
      en2 += x * tmp2 / tmp1;
      en3 += x * tmp3 / tmp1;
      llf += x * log(tmp1) - loggamma(x+1);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += log(t);
      en3 += log(t) * log(t);
      llf += log(lnorm_pdf(t, shape, scale));
    }
  }
  llf += log(omega) * en1;  // en1 is total number of faults
  en1 += omega * g10;  // g10 is the last time
  *total = en1;
  en2 += omega * g11;  // g11 is the last time
  en3 += omega * g12;  // g12 is the last time
  llf += - omega * (1.0 - g10);
  
  // M-step
  omega =  en1;
  scale = en2 / en1;
  shape = sqrt(en3 / en1 - scale*scale);

  pdiff[0] = omega - para[0];
  pdiff[1] = shape - para[1];
  pdiff[2] = scale - para[2];

  para[0] = omega;
  para[1] = shape;
  para[2] = scale;

  *retllf = llf;
}

