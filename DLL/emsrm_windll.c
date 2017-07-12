#include "emsrm_windll.h"
#include "emsrm/include/numlib.h"
#include "emsrm.h"

// em exp

int __stdcall srmdllExpEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_exp_emstep(dsize, time, num, type, 
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllExpMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	exp_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllExpRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	exp_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllExpInverseMVF
(double *value, double *time, int *npara, double *para) {
	exp_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllExpMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	exp_mttf(stime, h, ss, npara, para);
	return 0;
}


// em gamma

int __stdcall srmdllGammaEMstep
(int *dsize, double *time, double *num, int *type, int *divide,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_gamma_emstep(dsize, time, num, type, divide,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllGammaMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	gamma_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllGammaRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	gamma_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllGammaInverseMVF
(double *value, double *time, int *npara, double *para) {
	gamma_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllGammaMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	gamma_mttf(stime, h, ss, npara, para);
	return 0;
}

// em pareto

int __stdcall srmdllParetoEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_pareto_emstep(dsize, time, num, type,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllParetoMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	pareto_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllParetoRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	pareto_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllParetoInverseMVF
(double *value, double *time, int *npara, double *para) {
	pareto_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllParetoMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	pareto_mttf(stime, h, ss, npara, para);
	return 0;
}

// em tnorm

int __stdcall srmdllTNormEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_tnorm_emstep(dsize, time, num, type,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllTNormMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	tnorm_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllTNormRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	tnorm_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllTNormInverseMVF
(double *value, double *time, int *npara, double *para) {
	tnorm_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllTNormMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	tnorm_mttf(stime, h, ss, npara, para);
	return 0;
}

// em lnorm

int __stdcall srmdllLNormEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_lnorm_emstep(dsize, time, num, type,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllLNormMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	lnorm_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllLNormRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	lnorm_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllLNormInverseMVF
(double *value, double *time, int *npara, double *para) {
	lnorm_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllLNormMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	lnorm_mttf(stime, h, ss, npara, para);
	return 0;
}

// em tlogist

int __stdcall srmdllTLogistEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_tlogist_emstep(dsize, time, num, type,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllTLogistMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	tlogist_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllTLogistRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	tlogist_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllTLogistInverseMVF
(double *value, double *time, int *npara, double *para) {
	tlogist_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllTLogistMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	tlogist_mttf(stime, h, ss, npara, para);
	return 0;
}

// em llogist

int __stdcall srmdllLLogistEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_llogist_emstep(dsize, time, num, type,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllLLogistMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	llogist_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllLLogistRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	llogist_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllLLogistInverseMVF
(double *value, double *time, int *npara, double *para) {
	llogist_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllLLogistMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	llogist_mttf(stime, h, ss, npara, para);
	return 0;
}

// em txvmax

int __stdcall srmdllTXvMaxEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_txvmax_emstep(dsize, time, num, type,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllTXvMaxMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	txvmax_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllTXvMaxRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	txvmax_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllTXvMaxInverseMVF
(double *value, double *time, int *npara, double *para) {
	txvmax_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllTXvMaxMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	txvmax_mttf(stime, h, ss, npara, para);
	return 0;
}

// em lxvmax

int __stdcall srmdllLXvMaxEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_lxvmax_emstep(dsize, time, num, type,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllLXvMaxMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	lxvmax_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllLXvMaxRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	lxvmax_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllLXvMaxInverseMVF
(double *value, double *time, int *npara, double *para) {
	lxvmax_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllLXvMaxMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	lxvmax_mttf(stime, h, ss, npara, para);
	return 0;
}

// em txvmin

int __stdcall srmdllTXvMinEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_txvmin_emstep(dsize, time, num, type,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllTXvMinMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	txvmin_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllTXvMinRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	txvmin_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllTXvMinInverseMVF
(double *value, double *time, int *npara, double *para) {
	txvmin_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllTXvMinMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	txvmin_mttf(stime, h, ss, npara, para);
	return 0;
}

// em lxvmin

int __stdcall srmdllLXvMinEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total) {
  em_lxvmin_emstep(dsize, time, num, type,
		npara, para, pdiff, retllf, total);
  return 0;
}

int __stdcall srmdllLXvMinMVF
(int *dsize, double *time, double *mean, int *npara, double *para) {
	lxvmin_mvf_series(dsize, time, mean, npara, para);
	return 0;
}

int __stdcall srmdllLXvMinRate
(int *dsize, double *time, double *rate, int *npara, double *para) {
	lxvmin_rate_series(dsize, time, rate, npara, para);
	return 0;
}

int __stdcall srmdllLXvMinInverseMVF
(double *value, double *time, int *npara, double *para) {
	lxvmin_inverse_mvf(value, time, npara, para);
	return 0;
}

int __stdcall srmdllLXvMinMTTF
(double *stime, double *h, double *ss, int *npara, double *para) {
	lxvmin_mttf(stime, h, ss, npara, para);
	return 0;
}
