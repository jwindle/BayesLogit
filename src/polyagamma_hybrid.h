// Jesse Windle, 2019

#ifndef POLYAGAMMA_HYBRID
#define POLYAGAMMA_HYBRID

void rpg_gamma(double *x, double *n, double *z, int *num, int *trunc);
void rpg_devroye(double *x, int *n, double *z, int *num);
void rpg_alt(double *x, double *h, double *z, int* num);
void rpg_sp(double *x, double *h, double *z, int* num, int *iter);
void rpg_hybrid(double *x, double *h, double *z, int* num);

#endif
