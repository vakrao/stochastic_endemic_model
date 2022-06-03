#include "initparams.h"
#include "stoch_models.h"
#include<stdlib.h>
#include<omp.h>
#include<stdio.h>
#include<time.h>
#include<gsl/gsl_rng.h>  
#include<pthread.h>
//#include<omp.h>

#define NTHREADS 10

//runs the model, takes as input from the command line # of simulations to complete 
// ASSUMPTIONS: assumes that we want data for every possible vv 
// ASSUMPTIONS: look at variable initilizations in simconstants.h 
// Reccomended number of runs is 1000 
extern gsl_rng *r;
//main function run
int main(int argc, char* argv[]){
    pthread_t threads[NTHREADS];
    int thread_args[NTHREADS];
    int rc,i;
    if(argc != 4){
        printf("Missing parameters!");
        return 0;
    }
    //extern const int AGES;  
    //extern const double b;  
    //extern const double sigma; 
    //extern const double sigma_h; 
    //extern const double sigma_d; 
    //extern const double ve1;
    //extern const double variant_eff_2;
    //extern const double C1;
    //extern const double C2;
    //extern const int years; 
    //extern const double zeta;
    //extern const double ti_icu;
    //extern const int ft;
    //extern const double VC1; 
    //extern const double VC2;
    //extern const double variant_start;
    //extern const double R01; 
    //extern const double gamma_infec;
    //extern const double time_of_waning_natural;
    //extern const double time_of_immunity;
    //extern double variant_start_R02;

    const char *vv_title =  "../data/vv_vals.csv";
    struct ParameterSet p;
    initialize_named_params("one_strain.csv",&p);
    p.vv_values = initialize_unique_csv(11,vv_title,p.vv_values);
    time_t seconds;

    int run_number = atoi(argv[1]);
    int percent_number = atoi(argv[2]);
    int model_type = atoi(argv[3]);
    int vv_index = percent_number/10;
    // initializes psi, the vaccine coverage variable 
    // based on vaccinating more individuasl in the first 
    // year than any other year 
    //initialize gsl random environment
    char* dynamic_title = (char*) malloc(sizeof(char)*100);
    char* new_file = (char*)malloc(sizeof(char)*90);
    fprintf(stderr,"MODEL SIGMAQ1: %f \n",p.sigma_q1);
    fflush(stderr);
    fprintf(stderr,"MODEL SIGMAD1: %f \n",p.sigma_d1);
    fflush(stderr);
    //Vaccine configs, relates to all the different vaccine percentages
    int vax_percent = percent_number / 10;
    fprintf(stderr,"Vaccine value: %lf, Percentage: %f \n",p.vv_values[2],percent_number);
    fflush(stderr);
    p.ft = p.years * 365;
    fprintf(stderr,"FT: %d, YEARS: %d \n",p.ft,p.years);
    fflush(stderr);
    
    //Running each vaccine percentage for number given by sim_number
    for(int j = 0; j < run_number+1; j++){
       new_file = generate_names(vax_percent,j);
       stoch_model(p.vv_values[vv_index],j,new_file,p,model_type);
	   free(new_file);
     } 
	new_file = (char*)malloc(sizeof(char)*90);
    free(p.vv_values);
    free(dynamic_title);
    free(new_file);
    free(p.mu);
    seconds = time(NULL);
    printf("MODEL DONE IN: %ld \n",seconds);
}
