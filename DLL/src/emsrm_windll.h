
#include <windows.h>
#include <winsock.h>


// em exp

int __stdcall srmdllExpEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllExpMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllExpRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllExpInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllExpMTTF
(double *stime, double *h, double *ss, int *npara, double *para);


// em gamma

int __stdcall srmdllGammaEMstep
(int *dsize, double *time, double *num, int *type, int *divide,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllGammaMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllGammaRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllGammaInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllGammaMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

// em pareto

int __stdcall srmdllParetoEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllParetoMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllParetoRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllParetoInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllParetoMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

// em tnorm

int __stdcall srmdllTNormEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllTNormMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllTNormRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllTNormInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllTNormMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

// em lnorm

int __stdcall srmdllLNormEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllLNormMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllLNormRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllLNormInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllLNormMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

// em tlogist

int __stdcall srmdllTLogistEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllTLogistMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllTLogistRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllTLogistInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllTLogistMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

// em llogist

int __stdcall srmdllLLogistEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllLLogistMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllLLogistRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllLLogistInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllLLogistMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

// em txvmax

int __stdcall srmdllTXvMaxEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllTXvMaxMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllTXvMaxRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllTXvMaxInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllTXvMaxMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

// em lxvmax

int __stdcall srmdllLXvMaxEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllLXvMaxMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllLXvMaxRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllLXvMaxInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllLXvMaxMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

// em txvmin

int __stdcall srmdllTXvMinEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllTXvMinMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllTXvMinRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllTXvMinInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllTXvMinMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

// em lxvmin

int __stdcall srmdllLXvMinEMstep
(int *dsize, double *time, double *num, int *type,
 int *npara, double *para, double *pdiff, double *retllf, double *total);

int __stdcall srmdllLXvMinMVF
(int *dsize, double *time, double *mean, int *npara, double *para);

int __stdcall srmdllLXvMinRate
(int *dsize, double *time, double *rate, int *npara, double *para);

int __stdcall srmdllLXvMinInverseMVF
(double *value, double *time, int *npara, double *para);

int __stdcall srmdllLXvMinMTTF
(double *stime, double *h, double *ss, int *npara, double *para);

