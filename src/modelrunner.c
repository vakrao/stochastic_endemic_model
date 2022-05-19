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

int main(int argc, char* argv[]){
    pthread_t threads[NTHREADS];
    int thread_args[NTHREADS];
    int rc,i;
    if(argc != 2){
        printf("Need number of runs!");
        return 0;
    }
    extern const int AGES;  
    extern const double b;  
    extern const double sigma; 
    extern const double sigma_h; 
    extern const double sigma_d; 
    extern const double ve1;
    extern const double variant_eff_2;
    extern const double C1;
    extern const double C2;
    extern const int years; 
    extern const double zeta;
    extern const double ti_icu;
    extern const int ft;
    extern const double VC1; 
    extern const double VC2;
    extern const double variant_start;
    extern const double R01; 
    extern const double gamma_infec;
    extern const double time_of_waning_natural;
    extern const double time_of_immunity;
    extern double variant_start_R02;

    double* vv_values;
    const char *vv_title =  "../data/vv_vals.csv";
    vv_values = initialize_unique_csv(11,vv_title,vv_values);
    time_t seconds;

    int run_number = atoi(argv[1]);
    int percent_number = atoi(argv[2]);
    // initializes psi, the vaccine coverage variable 
    // based on vaccinating more individuasl in the first 
    // year than any other year 
    //initialize gsl random environment
    char* dynamic_title = (char*) malloc(sizeof(char)*100);
    char* new_file = (char*)malloc(sizeof(char)*90);

    //Vaccine configs, relates to all the different vaccine percentages
    int vax_percent = percent_number * 10;
    fprintf(stderr,"Vaccine value: %lf, Percentage: %d \n",vv_values[i],vax_percent);
    fflush(stderr);
    //Running each vaccine percentage for number given by sim_number
    int vv_index = percent_number/10;
    for(int j = 0; j < run_number+1; j++){
       new_file = generate_names(percent_number,j);
       stoch_model(vv_values[vv_index],j,new_file);
	   free(new_file);
     } 
	new_file = (char*)malloc(sizeof(char)*90);
    free(vv_values);
    free(dynamic_title);
    free(new_file);
    seconds = time(NULL);
    printf("MODEL DONE IN: %ld \n",seconds);
    
}
