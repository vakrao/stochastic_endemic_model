#ifndef INIT_PARAMS_H
#define INIT_PARAMS_H
void initialize_repeated_csv(int list_size,const char* filename, double* lst);
void read_contact_matrices(int list_size,const char* filename, double** lst);
double* initialize_unique_csv(int list_size,const char* filename, double* lst);
int initialize_named_params();
double* assign_hospitalization(float ti_icu, float* ICU_ratio,int ages);
char* generate_names(int stoch_number, int vv_value);
#endif