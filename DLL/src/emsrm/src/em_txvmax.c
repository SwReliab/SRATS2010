#include <math.h>
#ifndef _MSC_VER
#include "numlib.h"
#else
#include "../include/numlib.h"
#endif

/////////////////////////////////////////////
// EM method
/////////////////////////////////////////////

double txvmax_pdf(double t, double shape, double scale) {
  double y = exp(-(t-scale)/shape);
  double y0 = exp(scale/shape);
  return y*exp(-y)/shape/(1.0-exp(-y0));
}
	
double txvmax_cdf(double t, double shape, double scale) {
  double y = exp(-(t-scale)/shape);
  double y0 = exp(scale/shape);
  return 1.0-(1.0-exp(-y))/(1.0-exp(-y0));
}
	
double txvmax_mvf(double t, double omega, double shape, double scale) {
  return omega * txvmax_cdf(t, shape, scale);
}

double txvmax_cdf_v(double t, int npara, double *para) {
	return txvmax_cdf(t, para[1], para[2]);
}

void txvmax_inverse_mvf
(double *value, double *time, int *npara, double *para) {
  *time = findQuantile(*value/para[0], txvmax_cdf_v, *npara, para);
}

void txvmax_rate_series
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
    lambda[i] = omega * txvmax_pdf(t, shape, scale);
  }
}

double txvmax_creli(double x, double s, 
		   double omega, double shape, double scale,
		   double sval, double ffp) {
  double tmp;
  tmp = txvmax_mvf(x + s, omega, shape, scale) - sval;
  tmp = (exp(-tmp) - ffp) / (1.0 - ffp);
  return tmp;
}

void txvmax_mttf
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
  
  sval = txvmax_mvf(*stime, omega, shape, scale);
  ffp = exp(-(omega - sval));

  rev = sval - log(*r + (1.0-*r)*ffp);
  txvmax_inverse_mvf(&rev, &b, npara, para);
  b = b - *stime;
  
  // trapezoidal inte
  h = b - a;
  t = h * (txvmax_creli(a, *stime, omega, shape, scale, sval, ffp)
	   + txvmax_creli(b, *stime, omega, shape, scale, sval, ffp)) / 2.0;

  for (k = 1; k <= MTTF_MAXITE; k++) {
    n = 2 * n;
    h = h / 2.0;
    s = 0.0;
    for (i=1; i<=n-1; i+=2) {
      s += txvmax_creli(a + i * h, *stime, omega, shape, scale, sval, ffp);
    }
    tn = t/2.0 + h * s;
    if (fabs(tn - t) < INT_EPS*fabs(t)) {
      break;
    }
    t = tn;
  }
  *ss = tn;
}

void txvmax_mvf_series
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
    mean[i] = txvmax_mvf(t, omega, shape, scale);
  }
}

void em_txvmax_emstep
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
  y0 = exp(scale/shape);
  omega = omega / (1.0 - exp(-y0));

  en1 = 0.0;
  en2 = 0.0;
  en3 = 0.0;
  llf = 0.0;
  g00 = 1.0 - exp(-y0);
  g01 = exp(-scale/shape)*(1.0 - (1.0 + y0)*exp(-y0));
  g02 = 1.0 - exp(-y0) * (1.0 + y0 * log(y0));
  t = time[0];
  x = num[0];
  y = exp(-(t-scale)/shape);
  g10 = g00;
  g11 = g01;
  g12 = g02;
  g20 = 1.0 - exp(-y);
  g21 = exp(-scale/shape)*(1.0 - (1.0 + y)*exp(-y));
  g22 = 1.0 - exp(-y) * (1.0 + y * log(y));
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
    en2 += exp(-t/shape);
    en3 += (t-scale)/shape * (1.0-y);
    llf += log(txvmax_pdf(t, shape, scale));
  }
  for (j=1; j<*dsize; j++) {
    x = num[j];
    if (time[j] != 0.0) {
      t += time[j];
      y = exp(-(t-scale)/shape);
      g10 = g20;
      g11 = g21;
      g12 = g22;
      g20 = 1.0 - exp(-y);
      g21 = exp(-scale/shape)*(1.0 - (1.0 + y)*exp(-y));
      g22 = 1.0 - exp(-y) * (1.0 + y * log(y));
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
      en2 += exp(-t/shape);
      en3 += (t-scale)/shape * (1.0-y);
      llf += log(txvmax_pdf(t, shape, scale));
    }
  }
  llf += (log(omega) + log(g00))* en1;  // en1 is total number of faults
  *total = en1 + omega * g20;
  en1 += omega * (1.0 - g00 + g20);  // g00 is the first, g20 is the last
  en2 += omega * (exp(-scale/shape) - g01 
		  + g21);  // g01 is the first, g21 is the last
  en3 += omega * (1.0 - g02 + g22);  // g02 is the first, g22 is the last
  llf += - omega * (g00 - g20);
		
  // M-step
  //	shape = shape + Math.log(en3) - Math.log(en1);
  scale = -shape * log(en2/en1);
  shape = shape * (en3 / en1);
  y0 = exp(scale/shape);
  omega =  en1 * (1.0 - exp(-y0));
  // *total = omega;
		
  pdiff[0] = omega - para[0];
  pdiff[1] = shape - para[1];
  pdiff[2] = scale - para[2];

  para[0] = omega;
  para[1] = shape;
  para[2] = scale;

  *retllf = llf;
}
