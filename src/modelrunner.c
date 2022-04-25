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
    // initializes psi, the vaccine coverage variable 
    // based on vaccinating more individuasl in the first 
    // year than any other year 
    //initialize gsl random environment
    char* dynamic_title = (char*) malloc(sizeof(char)*100);
    char* dynamic_vv = (char*) malloc(sizeof(char)*100);
    char* csv_title = (char*)malloc(sizeof(char)*20);
    char* new_file = (char*)malloc(sizeof(char)*90);
    csv_title = ".csv";

    int vaccine_configs = 11;
    int starting_config = 0;
    //Vaccine configs, relates to all the different vaccine percentages
    	for(int i = starting_config; i < vaccine_configs; i++){
            fprintf(stderr,"Vaccine value: %lf",vv_values[i]);
            fflush(stderr);
    	    //Running each vaccine percentage for number given by sim_number
    	    //for(int j = 0; j < run_number+1; j++){
    	    //   new_file = generate_names(i,j);
    	    //   stoch_model(vv_values[i],j,new_file);
    	    //   fprintf(stderr,"FINISHED VV %d, RUN %d \n",vv_values[i],j);
    	    //   fflush(stderr);
	        //   free(new_file);
    	    //}
	        //ew_file = (char*)malloc(sizeof(char)*90);
    	}
   
    seconds = time(NULL);
    printf("MODEL DONE IN: %ld \n",seconds);


    
}
