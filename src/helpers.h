#ifndef _HELPERS_H
#define _HELPERS_H
double find_lambda(double q,int age, double sigma, double* I, double* VI, double** M, double* N,double* S);
double total(double* L);
double* ageing(double* L, int a);
double dynamic_vv(double* ABR, double* N, int VR, double ideal_vax_coverage);
//double poisson_draw(double mu,double max_value);
#endif
