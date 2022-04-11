#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simconstants.h"
#include<gsl/gsl_randist.h>  
#include<gsl/gsl_rng.h>  

//sigma = reduction of infectiousness
double find_lambda(double q,int age, double sigma, double* I, double* VI, double** M, double* N,double* S){
    double row_sum = 0;
    double infec_contact = 0; 
    double breakthrough_contact = 0;
   //         fprintf(stderr,"VI at 0: %lf \n",VI[0]);
   //         fflush(stderr);
     for(int i = 0; i < AGES; i++){
         infec_contact += (M[age][i]*(I[i]/N[i]));
     }
     for(int i = 0; i < AGES; i++){
         breakthrough_contact += (M[age][i]*(VI[i]/N[i]));
     }
     row_sum = q*(infec_contact + sigma*breakthrough_contact);
    return row_sum;
}

double* ageing(double* L, int a){
    double *new_L = (double*) malloc(AGES*sizeof(double));
    for(int i = 0; i < AGES; i++){
        if(i == 0){
            new_L[i] = 10;
        }
        else
        {
            new_L[i] = L[i-1];
        }
    }
    free(L);
    return new_L;
}

double total(double* L){
    double sum = 0;
    for(int i = 0; i<AGES; i++){
        sum += L[i];
    }
    return sum;
}

//double poisson_draw(double mu, double max_value){
//
//    double draw_value = max_value+1;
//    // gsl_rng *test; 
//    // const gsl_rng_type *T;
//    // gsl_rng_env_setup();
//    // T = gsl_rng_default;
//    // test = gsl_rng_alloc(T);
//    if(mu >= max_value){
//        return max_value;
//    }
//    while(draw_value > max_value ){
//        draw_value = gsl_ran_poisson(r,mu);
//    }
//    return draw_value;
//}




double dynamic_vv(double* ABR, double* G, int VS, double ideal_vax_coverage){
    float total_N = total(G);
    int counter = 1; 
    float vv_val= 0;
    float test_WVC = 0;
    for(int i =0; i < 13000; i++){
        float test_vv = i / 1000;
        for(int j = 0; j < AGES; j++){
            test_WVC += ABR[j]*total_N;
        }
        test_WVC = test_WVC*test_vv*VS;
        if(abs(test_vv - ideal_vax_coverage) < .0001){
            vv_val = test_vv;
            counter += 1;
            break;
        }
    }
    
}
