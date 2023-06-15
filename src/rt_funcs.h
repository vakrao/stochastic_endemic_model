#include "initparams.h"  
#ifndef RT_FUNCS_H
#define RT_FUNCS_H
float mod_rt_calc(double* S,double* I,double* R,double* V, double* N,double** M,double* mu, double* m,double q1,struct ParameterSet p);
float q_calc(double* S,double* I,double* R,double* V, double* N,double** M,double* mu, double* m,int R0,struct ParameterSet p);
#endif
