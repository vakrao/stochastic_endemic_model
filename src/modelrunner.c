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
// Reccomended number of runs is 100 
extern gsl_rng *r;
//main function run
int main(int argc, char* argv[]){
    pthread_t threads[NTHREADS];
    int thread_args[NTHREADS];
    int rc,i;
    if(argc != 5){
        printf("Missing parameters!");
        return 0;
    }

    const char *vv_title =  "../data/vv_vals.csv";
    struct ParameterSet p;
    char* file_name = (char*)malloc(sizeof(char)*2000);
    int run_number = atoi(argv[1]);
    int percent_number = atoi(argv[2]);
    int model_type = atoi(argv[3]);
    int vv_index = percent_number/10;
    file_name = argv[4];
    initialize_named_params(file_name,&p);
    p.vv_values = initialize_unique_csv(11,vv_title,p.vv_values);
    time_t seconds;

    // initializes psi, the vaccine coverage variable 
    // based on vaccinating more individuasl in the first 
    // year than any other year 
    //initialize gsl random environment
    char* dynamic_title = (char*) malloc(sizeof(char)*100);
    char* new_file = (char*)malloc(sizeof(char)*90);
    //Vaccine configs, relates to all the different vaccine percentages
    int vax_percent = percent_number / 10;
    p.ft = p.years * 365;




    const char *ifr_file =  "../data/ifr.csv";
    const char *vax_file =  "../data/dailyvax.csv";
    const char *m_file =  "../data/new_mort.csv";
    const char *n_file = "../data/us_pop.csv";
    const char *im_file = "../data/immigration_prop.csv";
    const char *overall_file = "../data/overall_contacts.csv";
    const char *icu_file = "../data/icu_ratio.csv";
    const char *school_file = "../data/school_contacts.csv";    
    p.m = (double*) malloc(p.AGES * sizeof(double));
    p.mu = (double*) malloc(p.AGES * sizeof(double));
    p.M = (double**) malloc(p.AGES * sizeof(double*));
    p.q1 = q_calc(double* S,double* I,double* R,double* V, double* N,double** M,double* mu, double* m,int R0,struct ParameterSet p)
    initialize_unique_csv(p.AGES,m_file,p.m);
    initialize_unique_csv(p.AGES,ifr_file,p.mu);
    read_contact_matrices(p.AGES, overall_file,p.M);

    //Running each vaccine percentage for number given by sim_number
    for(int j = 0; j < run_number+1; j++){

       new_file = generate_names(vax_percent,j);
       // stoch_model takes in a vv_value, iteration number, file_name to write to, parameters, and 
      // result format
      //model_type 0 -> ages and all data
      //model_type 1 -> ages and lessened data 
      //model_type 2-> age-agnostic data
       stoch_model(p.vv_values[vv_index],j,new_file,p,model_type);
	   free(new_file);
     } 
	new_file = (char*)malloc(sizeof(char)*90);
    free(p.vv_values);
    free(dynamic_title);
    free(new_file);
    free(p.mu);
    free(file_name);
}
