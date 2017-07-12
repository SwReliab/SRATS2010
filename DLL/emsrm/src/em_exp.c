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

double exp_pdf(double t, double rate) {
  return rate * exp(-rate * t);
}

double exp_cdf(double t, double rate) {
  return 1.0 - exp(-rate * t);
}

double exp_mvf(double t, double omega, double rate) {
  return omega * exp_cdf(t, rate);
}

double exp_quantile(double q, double rate) {
  return (-1.0/rate * log(1.0 - q));
}

void exp_inverse_mvf
(double *value, double *time, int *npara, double *para) {
  *time = exp_quantile(*value/para[0], para[1]);
}

void exp_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para) {
  int i;

  double t;
  double omega;
  double rate;

  if (*npara != 2) {
    return;
  }
  omega = para[0];
  rate = para[1];

  t = 0;
  for (i=0; i<*dsize; i++) {
    t += time[i];
    lambda[i] = omega * exp_pdf(t, rate);
  }
}

double exp_creli(double x, double s, double omega, double rate,
		 double sval, double ffp) {
  double tmp;
  tmp = exp_mvf(x + s, omega, rate) - sval;
  tmp = (exp(-tmp) - ffp) / (1.0 - ffp);
  return tmp;
}

void exp_mttf
(double *stime, double *r, double *ss, int *npara, double *para) {
  int i, k;
  int n = 1;
  double s, h, t, tn;
  double ffp, sval;
  double rev, b, a = 0.0;
  double omega, rate;

  if (*npara != 2) {
    *ss = -1.0;
    return;
  }
  omega = para[0];
  rate = para[1];
  
  sval = exp_mvf(*stime, omega, rate);
  ffp = exp(-(omega - sval));

  rev = sval - log(*r * (1.0 - ffp) + ffp);
  exp_inverse_mvf(&rev, &b, npara, para);
  b = b - *stime;
  
  // trapezoidal inte
  h = b - a;
  t = h * (exp_creli(a, *stime, omega, rate, sval, ffp)
	   + exp_creli(b, *stime, omega, rate, sval, ffp)) / 2.0;

  for (k = 1; k <= MTTF_MAXITE; k++) {
    n = 2 * n;
    h = h / 2.0;
    s = 0.0;
    for (i=1; i<=n-1; i+=2) {
      s += exp_creli(a + i * h, *stime, omega, rate, sval, ffp);
    }
    tn = t/2.0 + h * s;
    if (fabs(tn - t) < INT_EPS*fabs(t)) {
      break;
    }
    t = tn;
  }
  *ss = tn;
}

void exp_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para) {
  int i;

  double t;
  double omega;
  double rate;

  if (*npara != 2) {
    return;
  }
  omega = para[0];
  rate = para[1];

  t = 0;
  for (i=0; i<*dsize; i++) {
    t += time[i];
    mean[i] = exp_mvf(t, omega, rate);
  }
}

void em_exp_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {

  int j;
  double omega, rate;
  double t0;
  double x1, t1;
  double en1, en2;
  double tmp1, tmp2, llf;

  if (*npara != 2) {
    *retllf = 0.0;
    return;
  }
  omega = para[0];
  rate = para[1];

  // E-step
  t0 = 0.0;
  t1 = time[0];
  x1 = num[0];
  en1 = 0.0;
  en2 = 0.0;
  llf = 0.0;
  if (x1 != 0.0) {
    tmp1 = 1 - exp(-rate*t1);
    tmp2 = 1.0/rate - (t1 + 1/rate) * exp(-rate*t1);
    en1 = x1;
    en2 = x1 * tmp2 / tmp1;
    llf = x1 * log(tmp1) - loggamma(x1+1.0);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += t1;
    llf += log(rate) - rate*t1;
  }
  for (j=1; j<*dsize; j++) {
    if (time[j] != 0.0) {
      t0 = t1;
      t1 = t0 + time[j];
    }
    x1 = num[j];
    if (x1 != 0.0) {
      tmp1 = exp(-rate*t0) - exp(-rate*t1);
      tmp2 = (t0 + 1.0/rate) * exp(-rate*t0) 
	- (t1 + 1.0/rate) * exp(-rate*t1);
      en1 += x1;
      en2 += x1 * tmp2 / tmp1;
      llf += x1 * log(tmp1) - loggamma(x1+1.0);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += t1;
      llf += log(rate) - rate*t1;
    }
  }
  llf += log(omega) * en1;  // en1 is total number of faults
  en1 += omega * exp(-rate*t1);  // t1 is the last time
  *total = en1;
  en2 += omega * (t1 + 1.0/rate) * exp(-rate*t1); // t1 is the last time
  llf += - omega * (1.0 - exp(-rate*t1));
  
  // M-step
  omega = en1;
  rate = en1 / en2;

  pdiff[0] = omega - para[0];
  pdiff[1] = rate - para[1];
  para[0] = omega;
  para[1] = rate;

  *retllf = llf;
}

// test code
/* int main() { */
/*   int npara = 2; */
/*   double para[] = {12.5379740391877, 8.40477519171119e-02}; */
/*   double stime = 19.0; */
/*   double r = 0.001; */
  
/*   double ss; */

/*   exp_mttf(&stime, &r, &ss, &npara, para); */
/*   printf("%le\n", ss); */
/* } */
