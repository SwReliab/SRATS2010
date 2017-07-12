/* numlib */

#include <math.h>
#ifndef _MSC_VER
#include "numlib.h"
#else
#include "../include/numlib.h"
#endif

#define FACTMAX 20

#define FQ_EPS 1.0e-8
#define FQ_TOLZ 1.0e-15
#define FQ_MAXCNT 1024

static double PI = 3.14159265358979324; // pi
static double LOG_2PI = 1.83787706640934548;
static double LOG_PI = 1.14472988584940017; // log(pi)
static double LOG2PIOVER2 = 0.9189385332046727417803297364; // log(2pi) / 2

static int N = 8;
//static double B0 = 1;
//static double B1 = (-1.0 / 2.0);
static double B2 = ( 1.0 / 6.0);
static double B4 = (-1.0 / 30.0);
static double B6 = ( 1.0 / 42.0);
static double B8 = (-1.0 / 30.0);
static double B10 = ( 5.0 / 66.0);
static double B12 = (-691.0 / 2730.0);
static double B14 = ( 7.0 / 6.0);
static double B16 = (-3617.0 / 510.0);

static double nfact[FACTMAX+1] = {
		1.0,                        // 0
		1.0,                        // 1
		2.0,                        // 2
		6.0,                        // 3
		24.0,                       // 4
		120.0,                      // 5
		720.0,                      // 6
		5040.0,                     // 7
		40320.0,                    // 8
		362880.0,                   // 9
		3628800.0,                  // 10
		39916800.0,                 // 11
		479001600.0,                // 12
		6227020800.0,               // 13
		87178291200.0,              // 14
		1307674368000.0,            // 15
		20922789888000.0,           // 16
		355687428096000.0,          // 17
		6402373705728000.0,         // 18
		121645100408832000.0,       // 19
		2432902008176640000.0       // 20
};

static double lognfact[FACTMAX+1] = {
		0.0,
		0.0,
		0.6931471805599453,
		1.791759469228055,
		3.1780538303479458,
		4.787491742782046,
		6.579251212010101,
		8.525161361065415,
		10.60460290274525,
		12.801827480081469,
		15.104412573075516,
		17.502307845873887,
		19.987214495661885,
		22.552163853123425,
		25.19122118273868,
		27.89927138384089,
		30.671860106080672,
		33.50507345013689,
		36.39544520803305,
		39.339884187199495,
		42.335616460753485
};

double loggamma(double x) {
	double v, w;
	v = 1;
	while (x<N) { v *=x; x++; }
	w = 1 / (x * x);
	return ((((((((B16 / (16 * 15)) * w + (B14 / (14 * 13))) * w
			+ (B12 / (12 * 11))) * w + (B10 / (10 * 9))) * w
			+ (B8 / (8 * 7))) * w + (B6 / (6 * 5))) * w
			+ (B4 / (4 * 3))) * w + (B2 / (2 * 1))) / x
			+ 0.5 * LOG_2PI - log(v) - x + (x - 0.5) * log(x);
}

double gamma(double x) {
	if (x < 0) {
		return PI / (sin(PI * x) * exp(loggamma(1-x)));
	}
	return exp(loggamma(x));
}


double psi(double x) {
	double v, w;
	v = 0;
	while (x < N) { v += 1 / x; x++; }
	w = 1 / (x * x);
	v += ((((((((B16 / 16) * w + (B14 /14)) * w
			+ (B12 / 12)) * w + (B10 / 10)) * w
			+ (B8 / 8)) * w + (B6 / 6)) * w
			+ (B4 / 4)) * w + (B2 / 2)) * w + 0.5 / x;
	return log(x) - v;
}

double polygamma(int n, double x) {
	int k;
	double t, u, v, w;
	u = 1;
	for(k=1-n; k<0; k++) u *= k;
	v = 0;
	while (x<N) { v +=1 / pow(x, n+1); x++; }
	w = x * x;
	t = (((((((B16
			* (n + 15.0) * (n + 14) / (16 * 15 * w) + B14)
			* (n + 13.0) * (n + 12) / (14 * 13 * w) + B12)
			* (n + 11.0) * (n + 10) / (12 * 11 * w) + B10)
			* (n + 9.0) * (n + 8) / (10 * 9 * w) + B8)
			* (n + 7.0) * (n + 6) / (8 * 7 * w) + B6)
			* (n + 5.0) * (n + 4) / (6 * 5 * w) + B4)
			* (n + 3.0) * (n + 2) / (4 * 3 * w) + B2)
					* (n + 1.0) * n / (2 * 1 * w)
					+ 0.5 * n / x + 1;
	return u * (t / pow(x, n) + n * v);
}

double fact(int s) {
	if (s <= FACTMAX) {
		return nfact[s];
	} else {
		return exp(loggamma(1.0 + s));
	}
}

double logfact(int s) {
	if (s <= FACTMAX) {
		return lognfact[s];
	} else {
		return loggamma(1.0 + s);
	}
}

/* gamma dist */

double p_gamma(double a, double x, double loggamma_a) {
	int k;
	double result, term, previous;
	if (x >= 1+a) return 1 - q_gamma(a, x, loggamma_a);
	if (x == 0)   return 0;
	result = term = exp(a * log(x) - x - loggamma_a) / a;
	for (k=1; k<1000; k++) {
		term *= x / (a+k);
		previous = result;
		result += term;
		if (result == previous) return result;
	}
	return result;
}

double q_gamma(double a, double x, double loggamma_a) {
	int k;
	double result, w, temp, previous;
	double la = 1, lb = 1 + x - a;
	if (x < 1+a) return 1 - p_gamma(a, x, loggamma_a);
	w = exp(a * log(x) - x - loggamma_a);
	result = w/lb;
	for (k=2; k<1000; k++) {
		temp = ((k-1-a)*(lb-la)+(k+x)*lb)/k;
		la = lb;
		lb = temp;
		w *= (k-1-a)/k;
		temp = w/(la*lb);
		previous = result;
		result += temp;
		if (result == previous) return result;
	}
	return result;
}

/* Normal dist */

double p_normal(double x) {
	if (x >= 0)
		return 0.5 * (1 + p_gamma(0.5, 0.5 * x * x, LOG_PI / 2));
	else
		return 0.5 * q_gamma(0.5, 0.5 * x * x, LOG_PI / 2);
}

double q_normal(double x) {
	if (x >= 0)
		return 0.5 * q_gamma(0.5, 0.5 * x * x, LOG_PI / 2);
	else
		return 0.5 * (1 + p_gamma(0.5, 0.5 * x * x, LOG_PI / 2));
}

double d_normal(double x) {
	return 1.0/sqrt(2.0*PI)*exp(-x*x/2.0);
}

/* find quantile */

double findQuantile(double q,
		    double (*cdf)(double, int, double*),
		    int npara, double *para) {
  int cnt;
  double a, b, c;
  if (q >= 1.0 - FQ_TOLZ) return 0;
  a = 0.0;
  b = 1.0;
  cnt = 0;
  while ((*cdf)(b, npara, para) < q && cnt++ < FQ_MAXCNT) {
    a = b;
    b *= 2.0;
  }
  cnt = 0;
  c = (a+b)/2.0;
  while (fabs(a - b)/a > FQ_EPS && cnt++ < FQ_MAXCNT) {
    c = (a+b)/2.0;
    if ((*cdf)(c, npara, para) > q) {
      b = c;
    } else {
      a = c;
    }
  }
  return c;
}
