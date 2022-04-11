#include<gsl/gsl_rng.h>  
#ifndef _stoch_models_h
#define _stoch_models_h
double poisson_draw(gsl_rng *r,double mu,double max_value);
void stoch_model(double vv, int simulation_number,char* fileName); 
#endif
