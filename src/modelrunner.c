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
//struct ParameterSet{
//    const int AGES;  
//    const double b;  
//    const double sigma; 
//    const double sigma_h; 
//    const double sigma_d; 
//    const double ve1;
//    const double variant_eff_2;
//    const double C1;
//    const double C2;
//    const int years; 
//    const double zeta;
//    const double ti_icu;
//    const int ft;
//    const double VC1; 
//    const double VC2;
//    const double variant_start;
//    const double R01; 
//    const double gamma_infec;
//    const double time_of_waning_natural;
//    const double time_of_immunity;
//    double variant_start_R02;
//    double* vv_values;
//    const int AGES = 85; 
//    const double birth_rate = .012/365.0;
//    //reduction in infectiousness for multiple strains
//    const double sigma_red_i1= 0.5;
//    const double sigma_red_i2= 0.5;
//    //vaccine protection against infection 
//    const double sigma_i1 = 0.5;
//    const double sigma_i2 = 0.4;
//    //vaccine protection against hosptilaization 
//    const double sigma_h1 = 0.9; 
//    const double sigma_h2 = 0.9; 
//    //vaccine protection against death 
//    const double sigma_d1 = 0.95; 
//    const double sigma_d2 = 0.95; 
//    const double C1 = 0.5;
//    const double C2 = 0.9;
//    const int years = 20; 
//    //hospitalization waning recovery rate
//    const double zeta = 0.125;
//    const double ti_icu = 8;
//    const int ft = years*365;
//      const int ft = 10*365;
//    const double VC1 = 0.5; 
//    const double VC2 = 0.9;
//    //school and variant variables 
//    const double variant_start = -10;
//    const int school_spring = 150; 
//    const int school_break = 95; 
//    const int school_fall = 120;
//    const int vaccine_start = school_spring + school_break;
//    const int first_vax_seas_dur = 100;
//    const int perm_vax_seas_dur = 60;
//    const double R01 = 5; 
//    const double gamma = 0.1587;
//    // immunity values
//    const double time_of_waning_natural = 200;
//    const double time_of_immunity = 200;
//    const double variant_start_R02 = 0;
//};
//
int main(int argc, char* argv[]){
    pthread_t threads[NTHREADS];
    int thread_args[NTHREADS];
    int rc,i;
    struct ParameterSet params; 
    //if(argc != 4){
    //    printf("Need number of runs!");
    //    return 0;
    //}
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
    params.vv_values = initialize_unique_csv(11,vv_title,params.vv_values);
    time_t seconds;

    int run_number = atoi(argv[1]);
    int percent_number = atoi(argv[2]);
    int vv_index = percent_number/10;
    // initializes psi, the vaccine coverage variable 
    // based on vaccinating more individuasl in the first 
    // year than any other year 
    //initialize gsl random environment
    char* dynamic_title = (char*) malloc(sizeof(char)*100);
    char* new_file = (char*)malloc(sizeof(char)*90);
    fprintf(stderr,"PERCENT NUMBER: %d \n",percent_number);
    fflush(stderr);
    //Vaccine configs, relates to all the different vaccine percentages
    int vax_percent = percent_number / 10;
    fprintf(stderr,"Vaccine value: %lf, Percentage: %d \n",vv_values[vv_index],percent_number);
    fflush(stderr);
    //Running each vaccine percentage for number given by sim_number
    for(int j = 0; j < run_number+1; j++){
       new_file = generate_names(vax_percent,j);
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
