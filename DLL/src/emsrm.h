// emsrm header

#ifndef EM_SRM_H
#define EM_SRM_H

// exp

void em_exp_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void exp_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void exp_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void exp_inverse_mvf
(double *value, double *time, int *npara, double *para);

void exp_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// gamma

void em_gamma_emstep
(int *dsize, double *time, double *num, int *type, int *divide,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void gamma_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void gamma_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void gamma_inverse_mvf
(double *value, double *time, int *npara, double *para);

void gamma_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// pareto

void em_pareto_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void pareto_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void pareto_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void pareto_inverse_mvf
(double *value, double *time, int *npara, double *para);

void pareto_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// tnorm

void em_tnorm_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void tnorm_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void tnorm_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void tnorm_inverse_mvf
(double *value, double *time, int *npara, double *para);

void tnorm_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// lnorm

void em_lnorm_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void lnorm_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void lnorm_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void lnorm_inverse_mvf
(double *value, double *time, int *npara, double *para);

void lnorm_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// tlogist

void em_tlogist_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void tlogist_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void tlogist_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void tlogist_inverse_mvf
(double *value, double *time, int *npara, double *para);

void tlogist_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// llogist

void em_llogist_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void llogist_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void llogist_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void llogist_inverse_mvf
(double *value, double *time, int *npara, double *para);

void llogist_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// txvmax

void em_txvmax_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void txvmax_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void txvmax_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void txvmax_inverse_mvf
(double *value, double *time, int *npara, double *para);

void txvmax_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// lxvmax

void em_lxvmax_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void lxvmax_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void lxvmax_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void lxvmax_inverse_mvf
(double *value, double *time, int *npara, double *para);

void lxvmax_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// txvmin

void em_txvmin_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void txvmin_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void txvmin_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void txvmin_inverse_mvf
(double *value, double *time, int *npara, double *para);

void txvmin_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// lxvmin

void em_lxvmin_emstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

void lxvmin_mvf_series
(int *dsize, double *time, double *mean, int *npara, double *para);

void lxvmin_rate_series
(int *dsize, double *time, double *lambda, int *npara, double *para);

void lxvmin_inverse_mvf
(double *value, double *time, int *npara, double *para);

void lxvmin_mttf(double *stime, double *h, double *ss, int *npara, double *para);

// ks test

void isKSTest90
(int *dsize, double *time, double *num, int *type, double *mean, 
 int *result);

void isKSTest95
(int *dsize, double *time, double *num, int *type, double *mean, 
 int *result);

#endif
