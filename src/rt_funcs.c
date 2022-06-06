#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdbool.h>
#include <stdio.h>
#include <time.h>
#include <omp.h>
#include "helpers.h"
#include "initparams.h"
#include<gsl/gsl_randist.h>  
#include<gsl/gsl_rng.h>  
#include<gsl/gsl_matrix.h>  

int rt_calc(float* S,float* I,float* R,float* V,struct ParameterSet p){
    gsl_matrix *FL = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    gsl_matrix *G = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    bool change = false;
    int row_counter = 0; 
    int column_counter = 0;
    int adjusted_index = 0;

    // Setting values for flow OUT of infected subsystem.
    for(int i = 0; i < p.AGES; i++){
       if(change == true){
           float value = (p.gamma + m[adjusted_index] + mu[adjusted_index]*p.gamma);
           gsl_matrix_set(FL,i,adjusted_index,value);
           adjusted_index += 1;
           change = false;
       } 
       else{
           float value = (p.gamma + m[adjusted_index] + mu[adjusted_index]*p.gamma*p.sigma_1d);
           gsl_matrix_set(FL,i,adjusted_index,value);
           adjusted_index += 1;
           change = true;
       }
    }
    // creating gain matrix!
    for(int i = 0; i < p.AGES*2; i++){
        for(int j = 0; j < p.AGES*2;j++){
            float A = M[i][j]*(N[row_counter]/N[column_counter]); 
            float B = M[i][j]*(N[row_counter]/N[column_counter])*p.sigma_d1; 
            float C = V[row_counter]*(1-p.sigma_i1)*(M[row_counter][column_counter]/N[column_counter]); 
            float D = V[row_counter]*(1-p.sigma_i1)*p.sigma_q1*(M[row_counter][column_counter]/N[column_counter]); 
             if((i % 2 == 1) && (j % 2 == 1)){
                 gsl_matrix_set(G,i,j,A);
             } 
             if((i % 2 == 1) && (j % 2 == 0)){
                 gsl_matrix_set(G,i,j,B);
                 column_counter += 1;
             } 
             if((i % 2 == 0) && (j % 2 == 1)){
                 gsl_matrix_set(G,i,j,C);
             } 
             if((i % 2 == 0) && (j % 2 == 0)){
                 gsl_matrix_set(G,i,j,D);
                 column_counter += 1;
             } 
        }
        column_counter = 1;
        if(i % 2 == 0){
            row_counter += 1;
        }
    }

}



