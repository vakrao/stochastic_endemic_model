#include "initparams.h"
#include "helpers.h"
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
// Reccomended number of runs is 100 
extern gsl_rng *r;
//main function run
int main(int argc, char* argv[]){
    pthread_t threads[NTHREADS];
    int thread_args[NTHREADS];
    int rc,i;
    if(argc != 6){
        printf("Missing parameters!");
        return 0;
    }

    const char *vv_title =  "../params/vv_vals.csv";
    struct ParameterSet p;
    char* file_name = (char*)malloc(sizeof(char)*2000);
    char* folder_type = (char*)malloc(sizeof(char)*2000);
    int run_number = atoi(argv[1]);
    int percent_number = atoi(argv[2]);
    int model_type = atoi(argv[3]);
    int vv_index = percent_number/10;
    file_name = argv[4];
    folder_type = argv[5];
    initialize_named_params(file_name,&p);
    p.vv_values = initialize_unique_csv(11,vv_title,p.vv_values);
    time_t seconds;

    // initializes psi, the vaccine coverage variable 
    // based on vaccinating more individuals in the first 
    // year than any other year 
    // initialize gsl random environment
    char* dynamic_title = (char*) malloc(sizeof(char)*100);
    char* new_file = (char*)malloc(sizeof(char)*90);
    //Vaccine configs, relates to all the different vaccine percentages
    int vax_percent = percent_number / 10;
    p.ft = 20 * 365;


    const char *ifr_file =  "../params/ifr.csv";
    const char *vax_file =  "../params/dailyvax.csv";
    const char *m_file =  "../params/new_mort.csv";
    const char *n_file = "../params/us_pop.csv";
    const char *im_file = "../params/immigration_prop.csv";
    const char *overall_file = "../params/overall_contacts.csv";
    const char *icu_file = "../params/icu_ratio.csv";
    const char *school_file = "../params/school_contacts.csv";    
    double *raw_N0 = (double*) malloc(sizeof(double)*p.AGES);
    p.m = (double*) malloc(p.AGES * sizeof(double));
    p.mu = (double*) malloc(p.AGES * sizeof(double));
    p.M = (double**) malloc(p.AGES * sizeof(double*));
    raw_N0 = initialize_unique_csv(p.AGES,n_file,raw_N0);
    for(i=0;i<p.AGES;i++){
        p.N0 += raw_N0[i];
    }
    initialize_unique_csv(p.AGES,m_file,p.m);
    initialize_unique_csv(p.AGES,ifr_file,p.mu);
    read_contact_matrices(p.AGES, overall_file,p.M);
    fprintf(stderr,"printing rt");
    fflush(stderr);

    //Running each vaccine percentage for number given by sim_number
    for(int j = 0; j < run_number; j++){
       new_file = generate_names(vax_percent,j+1,folder_type);
      // stoch_model takes in a vv_value, iteration number, file_name to write to, parameters, and 
      // result format
      // model_type 0 -> ages and all data
      // model_type 1 -> ages and lessened data 
      // model_type 2-> age-agnostic data
//       stoch_model(p.vv_values[vv_index],j+1,"mod_birth.csv",p,model_type,percent_number);
    stoch_model(p.vv_values[vv_index],j+1,new_file,p,model_type,percent_number);
	   free(new_file);
	   new_file = (char*)malloc(sizeof(char)*90);
     } 
    free(p.vv_values);
    free(dynamic_title);
    free(new_file);
}
