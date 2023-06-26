#ifndef INIT_PARAMS_H
#define INIT_PARAMS_H
struct ParameterSet{
     int AGES;                                    
     int N0;
     //reduction in infectiousness for multiple strains
     double sigma_q1;
     double sigma_q2;
     //vaccine protection against infection 
     double sigma_i1;
     double sigma_i2;
     //vaccine protection against hosptilaization 
     double sigma_h1; 
     double sigma_h2; 
     //vaccine protection against death 
     double sigma_d1; 
     double sigma_d2; 
     //vaccine cross protection 
     double sigma_VC1; 
     double sigma_VC2; 
     double C1;
     double C2;
     int years; 
     //hospitalization waning recovery rate
     double zeta;
     double ti_icu;
     int ft;
     double VC1; 
     double VC2;
     //school and variant variables 
     double variant_start;
     int school_spring; 
     int school_break; 
     int school_fall;
     int vax_start;
     int first_vax_seas_dur;
     int perm_vax_seas_dur;
     double R01; 
     float q1;
     double R02;
     double gamma;
     // immunity values
     double time_of_waning_natural;
     double time_of_immunity;
     double* age_based_coverage;
     //  vv values + death rate
    double* vv_values;
    double* mu;
    double* m;
    double** M;
    double b; 
    double IFR_mod;
};
extern struct ParameterSet p;
void initialize_repeated_csv(int list_size,const char* filename, double* lst);
void read_contact_matrices(int list_size,const char* filename, double** lst);
double* initialize_unique_csv(int list_size,const char* filename, double* lst);
int initialize_named_params(const char* filename, struct ParameterSet *p);
double* assign_hospitalization(float ti_icu, float* ICU_ratio,int ages);
char* generate_names(int stoch_number, int vv_value, char* folder);
#endif
