#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include<gsl/gsl_rng.h>  
//Make .csv files for
//vv, ABR, TI_ICU, ICU_ratio, M, immigration_proprotion
//anything else...?


/*
Description: returns float-list based on the repeated csv format
2 columned .csv file, where left cell represents value and right cell represents number of times to repeat
@Params: 2 params, size, filename
    list_size :size of list to initialize
    filename: filename of csv to read
    Reads specific type of .csv which is formatted 
    with two columns, left side represents a value, right side represents number of repitions 

@Returns: 
    returns list of floats with parameter values
*/
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
extern const double gamma;
extern const double time_of_waning_natural;
extern const double time_of_immunity;
extern double variant_start_R02;
extern gsl_rng *r;

void initialize_repeated_csv(int list_size,const char* filename, double* lst){
    FILE* stream = fopen(filename,"r");
    char line[300];
    float value = 0;
    int reps = 0;
    int counter = 0;
    while(fgets(line, 1024, stream)){
        value = atof(strtok(line,","));
        reps =  atoi(strtok(NULL,","));
        for(int i = counter; i < reps+counter; i++){
            lst[i] = (float) value;
        }        
        counter = reps + counter;
    }
    fclose(stream);
}

void read_contact_matrices(int list_size,const char* filename, double** lst){
    FILE* stream = fopen(filename,"r");
    char line[100000];
    float value = 0;
    int reps = 0;
    int counter = 0;
    int i = 0;
    while(fgets(line, 1640, stream)){
        lst[i] = (double*) malloc(sizeof(double) * AGES);
        
        char* edptr;
        char* string_values = (strtok(line,","));
        int j = 0;
        while(string_values != NULL) {
            double ex;
            ex = strtod(string_values, &edptr);
            lst[i][j] = ex;
            j += 1; 
            string_values = strtok(NULL, ",");
        }
        i += 1;
    }
   fclose(stream);
}

double* initialize_unique_csv(int list_size,const char* filename,double *lst){
    FILE* stream = fopen(filename,"r");
    char line[300];
    float value = 0;
    int reps = 0;
    int counter = 0;
    double* new_lst = (double*) malloc(30*list_size*sizeof(double));
    for(int i =0; i < list_size*30;i++){
	new_lst[i] = 0.0; 
    }

    while(fgets(line, 1024, stream)){
	//char* new_val = (char*) malloc(sizeof(char)*1024);
        char* new_val = strtok(line,",");
        value = strtod(new_val,NULL);
        new_lst[counter] = (double) value;
        counter += 1;
    }
    fclose(stream);
    return new_lst;
}

int initialize_named_params(){
    FILE* stream = fopen("sim_params.csv","r");
    char line[300];
    float value = 0;
    int reps = 0;
    int counter = 0;
    while(fgets(line, 1024, stream)){
        value = atof(strtok(line,","));
        counter += 1;
    }
    fclose(stream);
    return 0;
}

/*
void ageing(double* M, int size){
    double* aged_M = (double*)malloc(85 * sizeof(double));
    aged_M[0] = 0;
    for(int i = 1; i < size; i++){
        aged_M[i] = M[i-1];
    }
    free(M);
    return;
}
*/

double* assign_hospitalization(float ti_icu, float* ICU_ratio,int ages){
    double* theta1 = (double *)malloc(ages*sizeof(double));
    for(int i = 0; i < ages; i++){
        theta1[i] = ti_icu*ICU_ratio[i];
    }
    return theta1;
}

char* generate_names(int vv_value, int stoch_number){
    char* new_string =  (char*) malloc(100*sizeof(char));
    sprintf(new_string,"noage_run_%d_%d0.csv",stoch_number,vv_value);
    return new_string;

}
