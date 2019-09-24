// Jesse Windle, 2019

#ifndef POLYAGAMMA_HYBRID
#define POLYAGAMMA_HYBRID


extern void rpg_gamma(double *x, double *n, double *z, int *num, int *trunc);
extern void rpg_devroye(double *x, int *n, double *z, int *num);
extern void rpg_alt(double *x, double *h, double *z, int* num);
extern void rpg_sp(double *x, double *h, double *z, int* num, int *iter);
extern void rpg_hybrid(double *x, double *h, double *z, int* num);



#endif
