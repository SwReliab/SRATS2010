// EM GAM

#include <math.h>
#ifndef _MSC_VER
#include "numlib.h"
#else
#include "../include/numlib.h"
#endif

/////////////////////////////////////////////
// Integral
/////////////////////////////////////////////

static double em_gamma_int_pi = 3.14159265358979324;
static double em_gamma_int_eps = 1.0e-8;

#define MAX_NN 300

int em_gamma_n;
double x[MAX_NN];
double w[MAX_NN];

double shape;
double rate;

double int_func(double x) {
  double y, tmp;
  y = rate * x;
  tmp = log(rate) + (shape-1.0) * log(y) - y - loggamma(shape);
  return log(x) * exp(tmp);
}

void em_gamma_makeW(int n) {
  int i, l, m;
  double p0, p1, p2;
  double q0, q1, q2;
  double tmp, dt;
  
  switch(n) {
  case 1:
    x[0] = 0.0;
    w[0] = 2.0;
    return;
  case 2:
    x[0] = sqrt(1.0/3.0);
    w[0] = 1.0;
    x[1] = -x[0];
    w[1] = w[0];
    return;
  case 3:
    x[0] = sqrt(0.6);
    w[0] = 5.0/9.0;
    x[1] = 0.0;
    w[1] = 8.0/9.0;
    x[2] = -x[0];
    w[2] = w[0];
    return;
  }
  
  m = n/2;
  for (i=0; i<m; i++) {
    tmp = cos((i+1.0-1.0/4.0)/(n+1.0/2.0)*em_gamma_int_pi);
    do {
      p1 = tmp;
      p2 = (3.0*tmp*tmp-1.0)/2.0;
      q1 = 1.0;
      q2 = 3.0*tmp;
      for (l=3; l<=n; l++) {
	p0 = p1;
	p1 = p2;
	p2 = ((2.0*l-1)*tmp*p1-(l-1)*p0)/l;
	q0 = q1;
	q1 = q2;
	q2 = ((2.0*l-1)*(tmp*q1+p1)-(l-1)*q0)/l;
      }
      dt = p2/q2;
      tmp = tmp - dt;
    } while(fabs(dt) > fabs(tmp)*em_gamma_int_eps);
    x[i] = tmp;
    w[i] = 2.0/(n*p1*q2);
  }
  if (n % 2 != 0) {
    x[n/2] = 0.0;
    tmp = (double) n;
    for (i=1; i<=m; i++)
      tmp = tmp*(0.5 - i)/i;
    w[n/2] = 2.0/(tmp*tmp);
  }
  for (i=0; i<m; i++) {
    x[n-1-i] = -x[i];
    w[n-1-i] = w[i];
  }
  return;
}

double em_gamma_int(double a, double b) {
  int i;
  double t1, t2;
  double x1, v, sum;
  
  t1 = (b - a)/2.0;
  t2 = (b + a)/2.0;
  sum = 0.0;
  
  for (i=0; i<em_gamma_n; i++) {
    x1 = t1*x[i] + t2;
    v = w[i]*int_func(x1);
    sum += v;
  }
  sum *= t1;
  return sum;
}

/////////////////////////////////////////////
// optimize
/////////////////////////////////////////////

double findshape_eps = 1.0e-8;
int findshape_maxcnt = 200;
	
double em_gamma_findshape(double init, double v) {
  int cnt;
  double a, b, c;
  a = init/2.0;
  b = init;
  cnt = 0;
  while (log(b) - psi(b) > v 
	 && cnt++ < findshape_maxcnt) {
    a = b;
    b *= 2.0;
  }
  while (log(a) - psi(a) <= v 
	 && cnt++ < findshape_maxcnt) {
    b = a;
    a = a/2.0;
  }
  cnt = 0;
  c = (a+b)/2.0;
  while (fabs(a - b)/a > findshape_eps
	 && cnt++ < findshape_maxcnt) {
    c = (a+b)/2.0;
    if (log(c) - psi(c) < v) {
      b = c;
    } else {
      a = c;
    }
  }
  return c;
}

/////////////////////////////////////////////
// EM method
/////////////////////////////////////////////

double gamma_pdf(double t, double shape, double rate) {
  double y;
  y = rate * t;
  return rate * pow(y, shape-1) * exp(-y) / exp(loggamma(shape));
}

double gamma_cdf(double t, double shape, double rate) {
  double y;
  y = rate * t;
  return p_gamma(shape, y, loggamma(shape));
}

double gamma_mvf(double t, double omega, double shape, double rate) {
  return omega * gamma_cdf(t, shape, rate);
}

double gamma_cdf_v(double t, int npara, double *para) {
	return gamma_cdf(t, para[1], para[2]);
}

void gamma_inverse_mvf
(double *value, double *time, int *npara, double *para) {
  *time = findQuantile(*value/para[0], gamma_cdf_v, *npara, para);
}

void gamma_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para) {
  int i;

  double t;

  double omega;
  double shape;
  double rate;

  if (*npara != 3) {
    return;
  }
  omega = para[0];
  shape = para[1];
  rate = para[2];

  t = 0;
  for (i=0; i<*dsize; i++) {
    t += time[i];
    lambda[i] = omega * gamma_pdf(t, shape, rate);
  }
}

double gamma_creli(double x, double s, 
		   double omega, double shape, double rate,
		   double sval, double ffp) {
  double tmp;
  tmp = gamma_mvf(x + s, omega, shape, rate) - sval;
  tmp = (exp(-tmp) - ffp) / (1.0 - ffp);
  return tmp;
}

void gamma_mttf
(double *stime, double *r, double *ss, int *npara, double *para) {
  int i, k;
  int n = 1;
  double s, h, t, tn;
  double ffp, sval;
  double rev, b, a = 0.0;
  double omega, shape, rate;

  if (*npara != 3) {
    *ss = 0.0;
    return;
  }
  omega = para[0];
  shape = para[1];
  rate = para[2];
  
  sval = gamma_mvf(*stime, omega, shape, rate);
  ffp = exp(-(omega - sval));

  rev = sval - log(*r + (1.0-*r)*ffp);
  gamma_inverse_mvf(&rev, &b, npara, para);
  b = b - *stime;
  
  // trapezoidal inte
  h = b - a;
  t = h * (gamma_creli(a, *stime, omega, shape, rate, sval, ffp)
	   + gamma_creli(b, *stime, omega, shape, rate, sval, ffp)) / 2.0;

  for (k = 1; k <= MTTF_MAXITE; k++) {
    n = 2 * n;
    h = h / 2.0;
    s = 0.0;
    for (i=1; i<=n-1; i+=2) {
      s += gamma_creli(a + i * h, *stime, omega, shape, rate, sval, ffp);
    }
    tn = t/2.0 + h * s;
    if (fabs(tn - t) < INT_EPS*fabs(t)) {
      break;
    }
    t = tn;
  }
  *ss = tn;
}

void gamma_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para) {
  int i;

  double t;
  double omega;
  double shape;
  double rate;

  if (*npara != 3) {
    return;
  }
  omega = para[0];
  shape = para[1];
  rate = para[2];

  t = 0;
  for (i=0; i<*dsize; i++) {
    t += time[i];
    mean[i] = gamma_mvf(t, omega, shape, rate);
  }
}

void em_gamma_emstep
(int *dsize, double *time, double *num, int *type, int *divide,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {

  int j;
  double t0;
  double x1, t1;
  double tmp1, tmp2, tmp3, tmp4, llf;
  double a0, a1;
  double gam10, gam11, gam20, gam21;
  double omega;

  double en1;
  double en2;
  double en3;

  if (*npara != 3) {
    *retllf = 0.0;
    return;
  }
  omega = para[0];
  shape = para[1];
  rate = para[2];

  em_gamma_n = *divide;
  em_gamma_makeW(em_gamma_n);

  // E-step
  a0 = loggamma(shape);
  a1 = loggamma(shape + 1.0);
  en1 = 0.0;
  en2 = 0.0;
  en3 = 0.0;
  llf = 0.0;

  t0 = 0.0;
  t1 = time[0];
  x1 = num[0];
  gam10 = 1.0;
  gam11 = 1.0;
  gam20 = q_gamma(shape, rate*t1, a0);
  gam21 = q_gamma(shape+1.0, rate*t1, a1);
  tmp3 = em_gamma_int(1.0e-10, t1);
  tmp4 = tmp3;
  if (x1 != 0.0) {
    tmp1 = gam10 - gam20;
    tmp2 = (shape/rate) * (gam11 - gam21);
    en1 += x1;
    en2 += x1 * tmp2 / tmp1;
    en3 += x1 * tmp3 / tmp1;
    llf += x1 * log(tmp1) - loggamma(x1+1.0);
  }
  if (type[0] == 1) {
    en1 += 1.0;
    en2 += t1;
    en3 += log(t1);
    llf += shape*log(rate) + (shape-1.0)*log(t1) 
      - rate*t1 - loggamma(shape);
  }
  for (j=1; j<*dsize; j++) {
    x1 = num[j];
    if (time[j] != 0.0) {
      t0 = t1;
      t1 = t0 + time[j];
      gam10 = gam20;
      gam11 = gam21;
      gam20 = q_gamma(shape, rate*t1, a0);
      gam21 = q_gamma(shape+1, rate*t1, a1);
      tmp3 = em_gamma_int(t0, t1);
      tmp4 += tmp3;
    }
    if (x1 != 0.0) {
      tmp1 = gam10 - gam20;
      tmp2 = (shape/rate) * (gam11 - gam21);
      en1 += x1;
      en2 += x1 * tmp2 / tmp1;
      en3 += x1 * tmp3 / tmp1;
      llf += x1 * log(tmp1) - loggamma(x1+1.0);
    }
    if (type[j] == 1) {
      en1 += 1.0;
      en2 += t1;
      en3 += log(t1);
      llf += shape*log(rate) + (shape-1.0)*log(t1) 
	- rate*t1 - loggamma(shape);
    }
  }
  llf += log(omega) * en1;  // en1 is total number of faults
  en1 += omega * gam20;  // gam20 is the last time
  *total = en1;
  en2 += omega * (shape/rate) * gam21;  // gam21 is the last time
  en3 += omega * (psi(shape) - log(rate) - tmp4);
  llf += - omega * (1.0 - gam20);
		
  // M-step
  omega =  en1;
  shape = em_gamma_findshape(shape, log(en2/en1)-en3/en1);
  rate = shape * en1 / en2;
		
  pdiff[0] = omega - para[0];
  pdiff[1] = shape - para[1];
  pdiff[2] = rate - para[2];

  para[0] = omega;
  para[1] = shape;
  para[2] = rate;

  *retllf = llf;
}
	
