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
#include<gsl/gsl_permutation.h>  
#include<gsl/gsl_matrix.h>  
#include<gsl/gsl_eigen.h>  
#include<gsl/gsl_blas.h>  
#include<gsl/gsl_linalg.h>  

float rt_calc(double* S,double* I,double* R,double* V, double* N,double** M,double* mu, double* m,double q1,struct ParameterSet p){
    gsl_matrix *FL = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    gsl_matrix *G = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    gsl_permutation *inv_perm = gsl_permutation_alloc(p.AGES*2);
    gsl_matrix *K = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    gsl_matrix_set_zero(FL);
    gsl_matrix_set_zero(G);
    gsl_matrix_set_zero(K);
    bool change = true;
    int row_counter = 0; 
    int column_counter = 0;
    int adjusted_index = 0;
    int s;
    double value = 0;

    // Setting values for flow OUT of infected subsystem.
    for(int i = 0; i < p.AGES*2; i++){
       if(change == true){
           value = (p.gamma + m[adjusted_index] + mu[adjusted_index]*p.gamma);
           gsl_matrix_set(FL,i,i,value);
           change = false;
       } 
       else{
           value = (p.gamma + m[adjusted_index] + mu[adjusted_index])*(1-p.sigma_d1)*(p.gamma);
           gsl_matrix_set(FL,i,i,value);
           adjusted_index += 1;
           change = true;
       }
    }
    // creating gain matrix!
    for(int i = 0; i < p.AGES*2; i++){
        for(int j = 0; j < p.AGES*2;j++){
            double A = M[row_counter][column_counter]*(N[row_counter]/N[column_counter]); 
            double B = M[row_counter][column_counter]*(N[row_counter]/N[column_counter])*p.sigma_d1; 
            double C = V[row_counter]*(1-p.sigma_i1)*(M[row_counter][column_counter]/N[column_counter]); 
            double D = V[row_counter]*(1-p.sigma_i1)*p.sigma_q1*(M[row_counter][column_counter]/N[column_counter]); 
             if((i % 2 == 0) && (j % 2 == 0)){
                 gsl_matrix_set(G,i,j,A);
             } 
             if((i % 2 == 0) && (j % 2 == 1)){
                 gsl_matrix_set(G,i,j,B);
                 column_counter += 1;
             } 
             if((i % 2 == 1) && (j % 2 == 0)){
                 gsl_matrix_set(G,i,j,C);
             } 
             if((i % 2 == 1) && (j % 2 == 1)){
                 gsl_matrix_set(G,i,j,D);
                 column_counter += 1;
             } 
        }
        column_counter = 0;
        if(i % 2 == 1){
            row_counter += 1;
        }
    }
    // take inverse now!
    gsl_linalg_LU_decomp(FL,inv_perm,&s);
    gsl_matrix *inv = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    gsl_linalg_LU_invert(FL,inv_perm,inv);
    gsl_permutation_free(inv_perm);
    gsl_matrix_scale(inv,1.0);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,inv,G,0.0,K);

    
    gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(p.AGES*2);
    gsl_vector *ev = gsl_vector_alloc(p.AGES*2);
    gsl_vector_set_zero(ev);
    gsl_eigen_symm(K, ev, w);
    double rt = gsl_vector_max(ev)*q1; 
    gsl_matrix_free(FL);
    gsl_matrix_free(G);
    gsl_matrix_free(K);
    gsl_matrix_free(inv);
   

    return rt;
}
float q_calc(double* S,double* I,double* R,double* V, double* N,double** M,double* mu, double* m,int R0,struct ParameterSet p){
    gsl_matrix *FL = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    gsl_matrix *G = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    gsl_permutation *inv_perm = gsl_permutation_alloc(p.AGES*2);
    gsl_matrix *K = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    gsl_matrix_set_zero(FL);
    gsl_matrix_set_zero(G);
    gsl_matrix_set_zero(K);
    bool change = true;
    int row_counter = 0; 
    int column_counter = 0;
    int adjusted_index = 0;
    int s;
    double value = 0;

    // Setting values for flow OUT of infected subsystem.
    for(int i = 0; i < p.AGES*2; i++){
       if(change == true){
           value = (p.gamma + m[adjusted_index] + mu[adjusted_index]*p.gamma);
           gsl_matrix_set(FL,i,i,value);
           change = false;
       } 
       else{
           value = (p.gamma + m[adjusted_index] + mu[adjusted_index])*(1-p.sigma_d1)*(p.gamma);
           gsl_matrix_set(FL,i,i,value);
           adjusted_index += 1;
           change = true;
       }
    }
    // creating gain matrix!
    for(int i = 0; i < p.AGES*2; i++){
        for(int j = 0; j < p.AGES*2;j++){
            double A = M[row_counter][column_counter]*(N[row_counter]/N[column_counter]); 
            double B = M[row_counter][column_counter]*(N[row_counter]/N[column_counter])*p.sigma_d1; 
            double C = V[row_counter]*(1-p.sigma_i1)*(M[row_counter][column_counter]/N[column_counter]); 
            double D = V[row_counter]*(1-p.sigma_i1)*p.sigma_q1*(M[row_counter][column_counter]/N[column_counter]); 
             if((i % 2 == 0) && (j % 2 == 0)){
                 gsl_matrix_set(G,i,j,A);
             } 
             if((i % 2 == 0) && (j % 2 == 1)){
                 gsl_matrix_set(G,i,j,B);
                 column_counter += 1;
             } 
             if((i % 2 == 1) && (j % 2 == 0)){
                 gsl_matrix_set(G,i,j,C);
             } 
             if((i % 2 == 1) && (j % 2 == 1)){
                 gsl_matrix_set(G,i,j,D);
                 column_counter += 1;
             } 
        }
        column_counter = 0;
        if(i % 2 == 1){
            row_counter += 1;
        }
    }
    // take inverse now!
    gsl_linalg_LU_decomp(FL,inv_perm,&s);
    gsl_matrix *inv = gsl_matrix_alloc(p.AGES*2,p.AGES*2);
    gsl_linalg_LU_invert(FL,inv_perm,inv);
    gsl_permutation_free(inv_perm);
    gsl_matrix_scale(inv,1.0);
    gsl_blas_dgemm(CblasNoTrans,CblasNoTrans,1.0,inv,G,0.0,K);

    
    gsl_eigen_symm_workspace *w = gsl_eigen_symm_alloc(p.AGES*2);
    gsl_vector *ev = gsl_vector_alloc(p.AGES*2);
    gsl_vector_set_zero(ev);
    gsl_eigen_symm(K, ev, w);
    gsl_eigen_symm_free(w);
    float q = R0/gsl_vector_max(ev); 
    return q;
}
