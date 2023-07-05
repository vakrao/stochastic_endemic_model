#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "simconstants.h"
#include<gsl/gsl_randist.h>  
#include<gsl/gsl_rng.h>  
#include "initparams.h"
#include "helpers.h"

//sigma = reduction of infectiousness
double find_lambda(double q,int age, double sigma, double* I, double* VI, double** M, double* N,double* S){
    double row_sum = 0;
    double contacts = 0;
    double pop = 0;
    double foi = 0;
    int z  = age/5;
    for(int j = 0; j < 17; j++){
        int start_age = (j*5);
        int end_age = start_age + 4;
        double ratio = 0;
        pop = 0;
        contacts = 0;
        // calculating population and contacts for 5 age block
        for(int k=start_age; k <= end_age; k++){
           contacts += ((sigma*VI[k]) + (I[k]));
           pop += N[k];
        }
        ratio = contacts/pop;
        // multiply by shorteend contact matrix
        row_sum += M[z][j]*(ratio);
     }
     foi = row_sum*q;
    return foi;
}


//Creating vaccination array 
//double* create_vaccine(int* psi,int vv,){
//    double age_based_coverage[AGES];
//    for(int i = 0; i < AGES; i++){
//	    age_based_coverage[i] = 0; 
//    }
//    //assigning vaccine percentages based on age
//    for(int i = 0; i < AGES; i++){
//        fprintf(stderr,"VV VALUE: %f \n",vv);
//        fflush(stderr);
//        if ( i <= 5 & i > 0 ){
//            age_based_coverage[i] = 0;
//        }
//        if ( i <= 12 & i > 5 ){
//            age_based_coverage[i] = (.646*vv)/365.0;
//        }
//        if ( i <= 17 & i > 12 ){
//            age_based_coverage[i] = (.553*vv)/365.0;
//        }
//        if ( i <= 49 & i > 17 ){
//            age_based_coverage[i] = (.384*vv)/365.0;
//        }
//        if ( i <= 64 & i > 49 ){
//            age_based_coverage[i] = (.506*vv)/365.0;
//        }
//        if ( i <= 85  & i > 64 ){
//            age_based_coverage[i] = (.698*vv)/365.0;
//        }
//    }
//    return age_based_coverage;
//}
//int* set_psi(int vaccine_start,int first_vsd,int perm_up,int perm_vp){
//    //setting psi value
//    for(int i = 0; i < years; i++){
//        for(int j = 0; j < 365; j++){
//            if( i == 0){
//                if(j < vaccine_start){
//                    psi[psi_counter] = 0; 
//                }
//                if(j >= vaccine_start && j < (vaccine_start + first_vsd)){
//                    psi[psi_counter] = 1; 
//                }
//                if(j >= (vaccine_start + first_vsd)){
//                    psi[psi_counter] = 0;
//                }
//            }
//            if( i > 0){
//                if(j < perm_up){
//                    psi[psi_counter] = 0; 
//                }
//                if(j >= perm_up && j < perm_vp){
//                    psi[psi_counter] = 1; 
//                }
//                if(j >= perm_vp){
//                    psi[psi_counter] = 0;
//                }
//            }
//            psi_counter += 1;
//        }
//    }
//
//
//} 

double* ageing(double* L,struct ParameterSet p){
    double *new_L = (double*) malloc(AGES*sizeof(double));
    new_L[0] = 0;
    for(int i = 1; i < p.AGES; i++){
        if(i == 84){
            new_L[i] = L[i-1] + L[i];
        }
        else{
            new_L[i] = L[i-1];
        }
//	fprintf(stderr,"new age value: %lf",new_L[i]);
//	fflush(stderr);
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



/*
Dynamic_vv calculates a new vv value based on the population at the time, the vaccine regime
the desired vaccine coverage percentage 
*/
double dynamic_vv(double* ABR, double* N, int VD, int ideal_vax_coverage){
    float total_N = total(N);
    int counter = 1; 
    double vv_val= 0;
    double test_WVC = 0;
    double base_val = 0;
    double vax_decimal = ideal_vax_coverage/100.0;
    double tN= 0;
    for(int i = 0; i < 85; i++){
        test_WVC += ABR[i]*N[i];
        tN += N[i];
    }
    base_val = test_WVC;
    double test_vv;
    double opt_val = 0;
    for(int i =1; i < 13000; i++){
        test_vv = (i*1.0)/1000.0;
        test_WVC = (test_vv*VD*test_WVC)/tN;
        double diff = (vax_decimal - test_WVC);
        if(diff < .00001 && test_WVC > 0){
            vv_val = test_vv;
            counter += 1;
            break;
        }
        test_WVC = base_val; 
    }
    return vv_val;
}
