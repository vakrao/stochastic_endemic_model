#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <stdio.h>
#include <time.h>
#include "helpers.h"
#include "initparams.h"
#include "helpers.h"
#include<gsl/gsl_randist.h>  
#include<gsl/gsl_rng.h>  

// Creae funcion ha akes int .csv file and creaes
// intpus
// Name hese ex files: simulaion
// Each simulaion will have a ceraint # of runs 
// Each simulaion will have .csv files associaed 
// vv is used o deerminte he vaccintaion level
// run number is used o name he files
// Need to add dynamic vax changing
gsl_rng *r;
void stoch_model(double vv, int run_number,char* fileName){
    FILE *fptr = fopen(fileName,"w");
    //fprintf(stderr,"STARTING MODEL! \n ");
    //fflush(stderr);

    const int AGES = 85; 
    const double birth_rate = .012/365.0;
    //reduction in infectiousness for multiple strains
    const double sigma_red_i1= 0.5;
    const double sigma_red_i2= 0.5;
    //vaccine protection against infection 
    const double sigma_i1 = 0.8;
    const double sigma_i2 = 0.8;
    //vaccine protection against hosptilaization 
    const double sigma_h1 = 0.9; 
    const double sigma_h2 = 0.9; 
    //vaccine protection against death 
    const double sigma_d1 = 0.95; 
    const double sigma_d2 = 0.95; 
    const double C1 = 0.5;
    const double C2 = 0.95;
    const int years = 20; 
    const double zeta = 0.45;
    const double ti_icu = 3;
    const int ft = years*365;
    const double VC1 = 0.8; 
    const double VC2 = 0.8;
    const double variant_start = 300;
    const double R01 = 5; 
    const double gamma = 0.1587;
    const double time_of_waning_natural = 365;
    const double time_of_immunity = 2*365;
    const double variant_start_R02 = 0;
    double age_based_coverage[AGES];
    srand(time(NULL));
    long value = rand()%100000;
    gsl_rng_env_setup();
    r = gsl_rng_alloc(gsl_rng_default);
    gsl_rng_set(r,value);
    int truncated_percent = trunc(vv);
    if(truncated_percent > 4 && truncated_percent < 11){
        truncated_percent = truncated_percent - 1;
    }
    else{
        truncated_percent = 10;
    }

    for(int i = 0; i < AGES; i++){
        if ( i <= 5 & i > 0 ){
            age_based_coverage[i] = 0;
        }
        if ( i <= 12 & i > 5 ){
            age_based_coverage[i] = (.646*vv)/365.0;
        }
        if ( i <= 17 & i > 12 ){
            age_based_coverage[i] = (.553*vv)/365.0;
        }
        if ( i <= 49 & i > 17 ){
            age_based_coverage[i] = (.384*vv)/365.0;
        }
        if ( i <= 64 & i > 49 ){
            age_based_coverage[i] = (.506*vv)/365.0;
        }
        if ( i <= 85  & i > 64 ){
            age_based_coverage[i] = (.698*vv)/365.0;
        }
    }
    double* m;
    double** cm_school = (double**) malloc(AGES * sizeof(double*));
    double** cm_overall = (double**) malloc(AGES * sizeof(double*));
    double** M = (double**) malloc(AGES * sizeof(double*));
    double *mu_i1;
    double *mu_i2=(double*) malloc(AGES * sizeof(double));
    double* VC = (double*) malloc(AGES * sizeof(double));
    double* theta = (double*) malloc(AGES * sizeof(double));
    double* ICU_raio =(double*) malloc(AGES * sizeof(double)); 
    int* psi = (int*) malloc(years*365*sizeof(int));
    int* school_time = (int*) malloc(years*sizeof(int));
    double* im_prop = (double*) malloc(AGES * sizeof(double));
    double* S = (double*) malloc(AGES * sizeof(double));
    double* I1 = (double*) malloc(AGES * sizeof(double));
    double* VI1 = (double*) malloc(AGES * sizeof(double));
    double* R1 = (double*) malloc(AGES * sizeof(double));
    double* VR1= (double*) malloc(AGES * sizeof(double));
    double* VR2= (double*) malloc(AGES * sizeof(double));
    double* H1= (double*) malloc(AGES * sizeof(double));
    double* V= (double*) malloc(AGES * sizeof(double));
    double* I2= (double*) malloc(AGES * sizeof(double));
    double* VI2= (double*) malloc(AGES * sizeof(double*));
    double* R2= (double*) malloc(AGES * sizeof(double*));
    double* H2= (double*) malloc(AGES * sizeof(double*));
    double* Xsi1= (double*) malloc(AGES * sizeof(double));
    double* Xsi2= (double*) malloc(AGES * sizeof(double));
    double* XIV1 = (double*) malloc(AGES * sizeof(double));
    double* XIV2 = (double*) malloc(AGES * sizeof(double));
    double* XI1 = (double*) malloc(AGES * sizeof(double));
    double* XI2 = (double*) malloc(AGES * sizeof(double));
    double* XD = (double*) malloc(AGES * sizeof(double));
    double* N;

    const char *ifr_ile =  "../data/ifr.csv";
    const char *vax_ile =  "../data/dailyvax.csv";
    const char *m_ile =  "../data/new_mort.csv";
    const char *n_ile = "../data/us_pop.csv";
    const char *im_ile = "../data/immigration_prop.csv";
    const char *overall_ile = "../data/overall_contacts.csv";
    const char *icu_file = "../data/icu_ratio.csv";
    int psi_counter = 0;
    for(int i = 0; i < years; i++){
        for(int j = 0; j < 365; j++){
            if( i == 0){
                if(j < 200){
                    psi[psi_counter] = 0; 
                }
                if(j >= 200 && j < 300 ){
                    psi[psi_counter] = 1; 
                }
                if(j >= 300){
                    psi[psi_counter] = 0;
                }
            }
            if( i > 0){
                if(j < 200){
                    psi[psi_counter] = 0; 
                }
                if(j >= 200 && j < 265 ){
                    psi[psi_counter] = 1; 
                }
                if(j >= 265){
                    psi[psi_counter] = 0;
                }
            }
            psi_counter += 1;
        }
    }

    // Two houghs o program model
    // Be specific abou names 
    // Be super-general abou names...? Use 3-d array 

    const float q1 = 0.0544;
//    const float q1 = 0;
    float R02 = 0;
    int t = 0;
    float q2  = 0;
    const float q2_value = 0;
    const double ime_vax_immuniy = years*2; 
    const double ime_wanintg_naural = years*2;

    double Xr1[AGES];
    double Y1[AGES];
    double Y2[AGES];
    double YV1[AGES];
    double YV2[AGES];
    double omega1[AGES];
    double omega2[AGES];
    double omega3[AGES];
    double omega4[AGES];
    double theta1[AGES];
    double theta2[AGES];
    double theta3[AGES];
    double theta4[AGES];
    double sv[AGES];
    double Xsh1[AGES];
    double Xsh2[AGES];
    double Xr1i2[AGES];
    double Xr2i1[AGES];
    double Xvvi1[AGES];
    double Xvvi2[AGES];
    double Yh1r1[AGES];
    double Yh2r2[AGES];
    double Xvr1vi2[AGES];
    double Xvr2vi1[AGES];
    double Di1[AGES];
    double DVi1[AGES];
    double DVi2[AGES];
    double Di2[AGES];
    double Dh1[AGES];
    double Dh2[AGES];
    double Vh1[AGES];
    double Vh2[AGES];
    double Vr1[AGES];
    double Vs[AGES];
    double vh1vr1[AGES];
    double vh2vr2[AGES];
    double p_vs[AGES];
    double r1v[AGES];
    double r2v[AGES];
    double IS[AGES];
    double II1[AGES];
    double IVI1[AGES];
    double IR1[AGES];
    double IH1[AGES];
    double IV[AGES];
    double II2[AGES];
    double IVI2[AGES];
    double IR2[AGES];
    double SD[AGES];
    double SB[AGES];
    double I1D[AGES];
    double I2D[AGES];
    double VI1D[AGES];
    double VI2D[AGES];
    double R2D[AGES];
    double R1D[AGES];
    double VR2D[AGES];
    double VR1D[AGES];
    double VD[AGES];
    double H2D[AGES];
    double H1D[AGES];
    double IH2[AGES];
    double D[AGES];
    double IVR1[AGES];
    double IVR2[AGES];
    double temp_transition = 0;
    double muIFR = 0;
    double ifr_i2_scale = 1;
    int new_yearly_imports = 100;

    // Initilizaing large datasets
    mu_i1 = initialize_unique_csv(AGES,ifr_ile,mu_i1);
    m = initialize_unique_csv(AGES,m_ile,m);
    N = initialize_unique_csv(AGES,n_ile,N);
    initialize_repeated_csv(AGES,vax_ile,VC);
    initialize_repeated_csv(AGES,im_ile,im_prop);
    initialize_repeated_csv(AGES,icu_file,ICU_raio);
    read_contact_matrices(AGES, overall_ile,cm_overall);
    int couner = 0;
    
    //fprintf(stderr,"entering AGE! \n");
    //fflush(stderr);
    for(int i = 0; i < AGES; i++){
//        theta[i] = (1/ti_icu)*ICU_raio[i];
        theta[i] = 0;
        mu_i1[i] = mu_i1[i] / 1000.0;
        mu_i2[i] = mu_i2[i] * ifr_i2_scale;
        m[i] = m[i] / 365.0;
        m[i] = m[i] * 5.0;
    }
    for(int c = 0; c < AGES; c++){
        M[c] = (double*) malloc(AGES*sizeof(double));
        S[c] = N[c];
        for(int j = 0; j < AGES; j++){
            M[c][j] = cm_overall[c][j];
        }
        VR2[c] = 1;
        I1[c] = 0;
        VI1[c] = 1; 
        R1[c] = 1;
        VR1[c] = 1;
        VR2[c] = 1;
        H1[c] = 1;
        V[c] = 0;
        I2[c] = 0; 
        VI2[c] = 0;
        R2[c] = 1; 
        H2[c] = 0; 
        Xsi1[c] = 1;
        Xsi2[c] = 0;
        Xr1[c] = 0;
        Y1[c] = 0;
        Y2[c] = 0;
        YV1[c] = 0;
        YV2[c] = 0;
        omega1[c] = 0;
        omega2[c] = 0;
        omega3[c] = 0;
        omega4[c] = 0;
        theta1[c] = 0;
        theta2[c] = 0;
        theta3[c] = 0;
        theta4[c] = 0;
        sv[c] = 0;
        Xsh1[c] = 0;
        Xsh2[c] = 0;
        Xr1i2[c] = 0;
        Xr2i1[c] = 0;
        Xvvi1[c] = 0;
        Xvvi2[c] = 0;
        Yh1r1[c] = 0;
        Yh2r2[c] = 0;
        Xvr1vi2[c] = 0;
        Xvr2vi1[c] = 0;
        Di1[c] = 0;
        DVi1[c] = 0;
        DVi2[c] = 0;
        Di2[c] = 0;
        Dh1[c] = 0;
        Dh2[c] = 0;
        Vh1[c] = 0;
        Vh2[c] = 0;
        Vr1[c] = 0;
        Vs[c] = 0;
        vh1vr1[c] = 0;
        vh2vr2[c] = 0;
        p_vs[c] = 0;
        r1v[c] = 0;
        r2v[c] = 0;
        IS[c] = 0;
        II1[c] = 0;
        IVI1[c] = 0;
        IR1[c] = 0;
        IH1[c] = 0;
        IV[c] = 0;
        II2[c] = 0; 
        IVI2[c] = 0;
        IR2[c] = 0;
        SD[c] = 0;
        SB[c] = 0;
        I1D[c] = 0;
        I2D[c] = 0;
        VI1D[c] = 0;
        VI2D[c] = 0;
        R2D[c] = 0;
        R1D[c] = 0;
        VR2D[c] = 0;
        VR1D[c] = 0;
        VD[c] = 0;
        H2D[c] = 0;
        H1D[c] = 0;
        IH2[c] = 0;
        D[c] = 0;
        IVR1[c] = 0;
        IVR2[c] = 0;
    }
   t = 0;
   fprintf(stderr,"STARTING TO ENTER LOOP!");
   fflush(stderr);
//    for(int c = 0; c < AGES; c++){
//        fprintf(stderr,"VI VALS: %f",VI1[c]);
//        fflush(stderr);
//    }
    double total_lambda = 0;
    int rand_number = 0;
    //Time loop starts here
    while(t < ft){
        if(t == variant_start){
            R02 = variant_start_R02;
            //To ADD: add new 100 imports ditributed across all ages
            for(int i = 0 ; i < new_yearly_imports; i++){
                rand_number = rand() % AGES;
                if(S[rand_number] > 0){
                    S[rand_number] = S[rand_number] -  1; 
                    I2[rand_number] = I2[rand_number] + 1; 
                }
            }
//            S[25] = S[25] - 1;
//            I2[25] = I2[25] + 1;
            q2 = q2_value;
        }
        else{
            if(t < variant_start && t == 0){
                for(int i = 0 ; i < new_yearly_imports; i++){
                    rand_number = rand() % AGES;
                    if(S[rand_number] > 0){
                        S[rand_number] = S[rand_number] - 1; 
                        I1[rand_number] = I2[rand_number] + 1; 
                    }
                }
               q2 = 0;  
            }
        }

        // Ageing Loop

        if(t % 365 == 0  & t != 0){
          S = ageing(S, AGES);
          I1 = ageing(I1, AGES);
          I2 = ageing(I2, AGES);
          VI1 = ageing(VI1, AGES);
          VI2 = ageing(VI2, AGES);
          V = ageing(V, AGES);
          R1= ageing(R1, AGES);
          R2= ageing(R2, AGES);
          H1= ageing(H1, AGES);
          H2= ageing(H2, AGES);
          VR1= ageing(VR1, AGES);
          VR2= ageing(VR2, AGES);
          S[0] = 1;
         if(t < variant_start){
            for(int i = 0 ; i < new_yearly_imports; i++){
                rand_number = rand() % AGES;
                if(S[rand_number] > 0){
                    S[rand_number] -= 1; 
                    I1[rand_number] += 1; 
                }
            }
         }
         else{
            for(int i = 0 ; i < new_yearly_imports; i++){
                rand_number = rand() % AGES;
                if(S[rand_number] > 0){
                    S[rand_number] -= 1; 
                    I2[rand_number] += 1; 
                }
            }
         }
         for(int k=0; k < AGES; k++){
             N[k] = S[k] + I1[k] + I2[k]+ VI1[k]+ VI2[k]+ V[k]+ R1[k]+ R2[k]+ H1[k]+ H2[k]+ VR1[k]+ VR2[k];
         }

        }
        else{
            if(t  == 0){
                S[25] = S[25] - 1;
                I1[25] = I1[25] + 1; 
            }
        }

        // Stochasic Age-Transition Loop
        for(int i=0; i < AGES; i++){
            double lambda1 = find_lambda(q1,i,sigma_red_i1,I1,VI1,M,N,S);
            total_lambda += lambda1;
            double lambda2 = find_lambda(q2,i,sigma_red_i2,I2,VI2,M,N,S);
            double lambda_mean = (S[i] * lambda1);
            Xsi1[i] =  poisson_draw(lambda_mean,S[i]);
            if(S[i] - Xsi1[i] >= 1){
                temp_transition = (S[i] - Xsi1[i]);
                Xsi2[i] = poisson_draw(temp_transition*lambda2,temp_transition);
            }
            else{
                Xsi2[i] = 0;
            }

            if(S[i] - Xsi1[i] - Xsi2[i] >= 1){
                temp_transition = (S[i] - Xsi1[i] - Xsi2[i]);
                sv[i] = poisson_draw(temp_transition*psi[t]*age_based_coverage[i],temp_transition);
            }
            else{
                sv[i] = 0;
            }

            if(S[i] - Xsi1[i] - Xsi2[i] - sv[i] >= 1){
                temp_transition = (S[i] - Xsi1[i] - Xsi2[i] - sv[i]);
                SD[i] = poisson_draw(temp_transition * m[i],temp_transition);
           }
            else{
                SD[i] = 0;
           }
           // V comparment exit logic
           Vs[i] = poisson_draw(V[i]*(1/time_of_immunity),V[i]);

            // sigma_i1 = sigma_i1
            // ve2 = sigma_i2
            if(V[i] - Vs[i] >= 1){
                temp_transition = V[i] - Vs[i];
                Xvvi1[i] = poisson_draw(temp_transition * lambda1 * (1-sigma_i1),temp_transition);
            }
            else{
                Xvvi1[i] = 0;
            }
            // fprintf(stderr,"Xvvi1: %lf \n",Xvvi1[i]);
            // fflush(stderr);

            if(V[i] - Vs[i] - Xvvi1[i] >= 1){
                temp_transition = V[i] - Vs[i] - Xvvi1[i];
                Xvvi2[i] = poisson_draw(temp_transition * lambda2 * (1-sigma_i2),temp_transition);
            }
            else{
                Xvvi2[i] = 0;
            }

            if(V[i] - Vs[i] - Xvvi1[i] - Xvvi2[i]>= 1){
                temp_transition = V[i] - Vs[i] - Xvvi1[i] - Xvvi2[i];
                DVi2[i] = poisson_draw(temp_transition * mu_i1[i],temp_transition);
            }
            else{
                DVi2[i] = 0;
            }

            if(V[i] - Vs[i] - Xvvi1[i] - Xvvi2[i] - DVi2[i] >= 1){
                temp_transition = V[i] - Vs[i] - Xvvi1[i] - Xvvi2[i] - DVi2[i];
                VI2D[i] = poisson_draw(temp_transition * mu_i2[i],temp_transition);
            }
            else{
                VI2D[i] = 0;
            }

            if(I1[i] >= 1){
                Y1[i] = poisson_draw(I1[i]*gamma,I1[i]);
            }
            else{
                Y1[i] = 0;
            }

            if(I1[i] - Y1[i] >= 1){
                temp_transition = (I1[i] - Y1[i]);
                theta1[i] = poisson_draw(temp_transition*theta[i],temp_transition);
            }
            else{
                theta1[i] = 0;
            }

            if(I1[i] - Y1[i] - theta1[i] >= 1){
                temp_transition = (I1[i] - Y1[i] - theta1[i]);
                Di1[i] = poisson_draw(temp_transition*mu_i1[i]*gamma,temp_transition);
            }
            else{
                Di1[i] = 0;
            }

            if(I1[i] - Y1[i] - theta1[i] - Di1[i] >= 1){
                temp_transition = (I1[i] - Y1[i] - theta1[i] - Di1[i]);
                I1D[i] = poisson_draw(temp_transition*m[i],temp_transition);
            }
            else{
                I1D[i] = 0;
            }

            // I2 exit logic, recovered
            if(I2[i] >= 1){
                Y2[i] = poisson_draw(I2[i]*gamma,I2[i]);
            }
            else{
                Y2[i] = 0;
            }
            // fprintf(stderr,"Y2[i]: %lf \n",Y2[i]);
            // fflush(stderr);

            // I2 exit logic, hospitalization
            if(I2[i] - Y2[i] >= 1){
                temp_transition = (I2[i] - Y2[i] );
                theta2[i] = poisson_draw(temp_transition*theta[i],temp_transition);
            }
            else{
                theta2[i] = 0;
            }

            if(I2[i] - Y2[i] - theta2[i] >= 1){
                temp_transition = (I2[i] - Y2[i] - theta2[i] );
                Di2[i] = poisson_draw(temp_transition*mu_i2[i]*gamma,temp_transition);
            }
            else{
                Di2[i] = 0;
            }
            // fprintf(stderr,"Di2[i]: %lf \n",Di2[i]);
            // fflush(stderr);

            if(I2[i] - Y2[i] - theta2[i] - Di2[i] >= 1 ){
                temp_transition = (I2[i] - Y2[i] - theta2[i] );
                I2D[i] = poisson_draw(temp_transition*m[i],temp_transition);
            }
            else{
                I2D[i] = 0;
            }

            // VI1 exit logic

            YV1[i] = poisson_draw(VI1[i]*gamma,VI1[i]);

            if(VI1[i] - YV1[i] >= 1){
                temp_transition = (VI1[i] - YV1[i] );
                theta3[i] = poisson_draw(temp_transition*theta[i]*(1-sigma_h1),temp_transition);
            }
            else{
                theta3[i] = 0;
            }

            if(VI1[i] - YV1[i] - theta3[i] >= 1){
                temp_transition = (VI1[i] - YV1[i] - theta3[i] );
                DVi1[i] = poisson_draw(temp_transition*theta[i]*(1-sigma_d1)*gamma,temp_transition);
            }
            else{
                DVi1[i] = 0;
            }

            if(VI1[i] - YV1[i] - theta3[i] - DVi1[i] >= 1){
                temp_transition = (VI1[i] - YV1[i] - theta3[i] - DVi1[i]);
                VI1D[i] = poisson_draw(temp_transition*m[i],temp_transition);
            }
            else{
                VI1D[i] = 0;
            }

            YV2[i] = poisson_draw(VI2[i]*gamma,VI2[i]);

            if(VI2[i] - YV2[i] >= 1){
                temp_transition = (VI2[i] - YV2[i]);
                theta4[i] = poisson_draw(temp_transition*theta[i]*sigma_h2,temp_transition);
            }
            else{
                theta4[i] = 0;
            }
            if(VI2[i] - YV2[i] - theta4[i] >= 1){
                temp_transition = (VI2[i] - YV2[i] - theta4[i]);
                DVi2[i] = poisson_draw(temp_transition*mu_i1[i]*gamma*sigma_d2,temp_transition);
            }
            else{
                DVi2[i] = 0;
            }
            if(VI2[i] - YV2[i] - theta4[i] - DVi2[i] >= 1){
                temp_transition = (VI2[i] - YV2[i] - theta4[i]);
                VI2D[i] = poisson_draw(temp_transition*m[i],temp_transition);
            }
            else{
                VI2D[i] = 0;
            }

            // Recovered Exit Logic to Susceptible
            omega1[i] = poisson_draw(R1[i]*(1/time_of_waning_natural),R1[i]);

            if(R1[i] - omega1[i]>= 1){
                temp_transition = R1[i] - omega1[i];
                r1v[i] = poisson_draw(temp_transition*psi[t]*age_based_coverage[i],temp_transition);
            }
            else{
                r1v[i] = 0;
            }

            if(R1[i] - omega1[i] - r1v[i]>= 1){
                temp_transition = R1[i] - omega1[i] -r1v[i];
                Xr1i2[i] = poisson_draw(temp_transition*(1-C2)*lambda2,temp_transition);
            }
            else{
                Xr1i2[i] = 0;
            }

            if(R1[i] - omega1[i] - r1v[i] - Xr1i2[i] >= 1){
                temp_transition = R1[i] - omega1[i] -r1v[i] - Xr1i2[i];
                R1D[i] = poisson_draw(temp_transition*m[i],temp_transition);
            }
            else{
                R1D[i] = 0;
            }

            //R2 exit logic
            omega2[i] = poisson_draw(R2[i]*(1/time_of_waning_natural),R2[i]);

            if(R2[i] - omega2[i]>= 1){
                temp_transition = R2[i] - omega2[i];
                r2v[i] = poisson_draw(temp_transition*psi[t]*age_based_coverage[i],temp_transition);
            }
            else{
                r2v[i] = 0;
            }

            if(R2[i] - omega2[i] - r2v[i]>= 1){
                temp_transition = R2[i] - omega2[i] -r2v[i];
                Xr2i1[i] = poisson_draw(temp_transition*(1-C2)*lambda1,temp_transition);
            }
            else{
                Xr1i2[i] = 0;
            }

            if(R2[i] - omega2[i] - r2v[i] - Xr2i1[i] >= 1){
                temp_transition = R2[i] - omega2[i] -r2v[i] - Xr2i1[i];
                R2D[i] = poisson_draw(temp_transition*m[i],temp_transition);
            }
            else{
                R2D[i] = 0;
            }

           Xvr1vi2[i] = poisson_draw(VR1[i]*(1-VC1)*lambda2*sigma_red_i1,VR1[i]); 

           if(VR1[i] - Xvr1vi2[i] >= 1){
               temp_transition = VR1[i] - Xvr1vi2[i];
               omega3[i] = poisson_draw(temp_transition*(1/time_of_immunity),temp_transition);
           }
           else{
               omega3[i] = 0;
           }

           if(VR1[i] - Xvr1vi2[i] - omega3[i] >= 1){
               temp_transition = VR1[i] - Xvr1vi2[i] - omega3[i];
               VR1D[i] = poisson_draw(temp_transition*m[i],temp_transition);
           }
           else{
               VR1D[i] = 0;
           }

           Xvr2vi1[i] = poisson_draw(VR2[i]*(1-VC2)*lambda1*sigma_red_i2,VR2[i]); 

           if(VR2[i] - Xvr2vi1[i] >= 1){
               temp_transition = VR2[i] - Xvr2vi1[i];
               omega4[i] = poisson_draw(temp_transition*(1/time_of_immunity),temp_transition);
           }
           else{
               omega4[i] = 0;
           }

           if(VR2[i] - Xvr2vi1[i] - omega4[i] >= 1){
               temp_transition = VR2[i] - Xvr2vi1[i] - omega4[i];
               VR2D[i] = poisson_draw(temp_transition*m[i],temp_transition);
           }
           else{
               VR2D[i] = 0;
           }

           if(H1[i] >= 1){
               Yh1r1[i] = poisson_draw(H1[i]*zeta,H1[i]);
           }
           else{
               Yh1r1[i] = 0;
           }

           if(H1[i] - Yh1r1[i] >= 1){
               temp_transition = H1[i]-Yh1r1[i];
               H1D[i] = poisson_draw(temp_transition*m[i],temp_transition);
           }
           else{
               H1D[i] = 0;
           }

           if(H2[i] >= 1){
               temp_transition = H2[i]-Yh2r2[i];
               Yh2r2[i] = poisson_draw(temp_transition*zeta,temp_transition);
           }
           else{
               Yh2r2[i] = 0;
           }

           if(H2[i] - Yh2r2[i] >= 1){
               temp_transition = H2[i]-Yh2r2[i];
               H2D[i] = poisson_draw(temp_transition*m[i],temp_transition);
           }
           else{
               H2D[i] = 0;
           }

            IS[i] =  poisson_draw((S[i] / N[i]) * im_prop[i],S[i]);
            IV[i] =  poisson_draw((V[i] / N[i]) * im_prop[i],V[i]);
            II1[i] = poisson_draw((I1[i] / N[i]) * im_prop[i],I1[i]);
            IH1[i] = poisson_draw(((H1[i] / N[i])) * im_prop[i],H1[i]);
            IR1[i] = poisson_draw((R1[i] / N[i]) * im_prop[i],R1[i]);
            IVI1[i]= poisson_draw((VI1[i] / N[i]) * im_prop[i],VI1[i]);
            IH2[i] = poisson_draw((H2[i] / N[i]) * im_prop[i],H2[i]);
            IR2[i] = poisson_draw((R2[i] / N[i]) * im_prop[i],R2[i]);
            IVI2[i]= poisson_draw((VI2[i] / N[i]) * im_prop[i],VI2[i]);
            IVR2[i]= poisson_draw((VR2[i] / N[i]) * im_prop[i],VR2[i]);

            if(i == 0){
                float totalN = total(N);
                SB[i] = poisson_draw(totalN*birth_rate,totalN);
            }
        }

        if(t == 0){
            fprintf(fptr,"t,vv,age,value,vartype\n");
        }

//        if(t == 0){
//            fprintf(fptr,"t,vv,value,vartype\n");
//        }


        for(int i = 0; i < AGES; i++){
            S[i] = S[i] + Vs[i] + omega1[i] + omega2[i] + omega3[i] + omega4[i] - Xsi1[i] - Xsi2[i] - sv[i] - SD[i] + SB[i] + IS[i];
            I1[i] = I1[i] + Xsi1[i] + Xr2i1[i] - Y1[i] - theta1[i] - Di1[i] - I1D[i] + II1[i];
            I2[i]  = I2[i] + Xsi2[i] + Xr1i2[i] - Y2[i] - theta2[i] - Di2[i] - I2D[i] + II2[i];
            R1[i]  = R1[i] + Y1[i] + Yh1r1[i] - r1v[i] - omega1[i] - Xr1i2[i] - R1D[i] + IR1[i];
            R2[i]  = R2[i] + Y2[i] + Yh2r2[i] - r2v[i] - omega2[i] - Xr2i1[i] - R2D[i] + IR2[i];
            VI1[i] = VI1[i] + Xvvi1[i] + Xvr2vi1[i] - theta3[i] - YV1[i] - DVi1[i] - VI1D[i] + IVI1[i];
            VI2[i] = VI2[i] + Xvvi2[i] + Xvr1vi2[i] - theta4[i] - YV2[i] - DVi2[i] - VI2D[i] + IVI2[i];
            VR1[i] = VR1[i] + YV1[i] - omega3[i] - Xvr1vi2[i] - VR1D[i] + IVR1[i];
            VR2[i] = VR2[i] + YV2[i] - omega4[i] - Xvr2vi1[i] - VR2D[i] + IVR2[i];
            V[i]   = V[i] + r1v[i] + r2v[i] + sv[i] - Xvvi1[i]  - Xvvi2[i] - Vs[i]  - VD[i] + IV[i];
            H1[i]  = H1[i] + theta1[i] + theta3[i] - Yh1r1[i] - H1D[i] + IH1[i];
            H2[i]  = H2[i] + theta2[i] + theta4[i] - Yh2r2[i] - H2D[i] + IH2[i];
            D[i]   = SD[i] + I1D[i] + I2D[i] + R1D[i] + R2D[i] + VD[i] + VI1D[i] + VI2D[i] + VR1D[i] + VR2D[i] + H1D[i] + H2D[i]; 
            N[i]   = S[i] + I1[i] + I2[i] + R1[i] + R2[i] + VI1[i] + VI2[i] + VR1[i] + VR2[i] + V[i] + H1[i] + H2[i];
            XIV1[i] = Xvvi1[i] + Xvr2vi1[i];
            XIV2[i] = Xvvi2[i] + Xvr1vi2[i];
            XI1[i] = Xsi1[i] + Xr2i1[i]; 
            XI2[i] = Xsi2[i] + Xr1i2[i];
            XD[i] = D[i] + Di1[i] + Di2[i] + DVi1[i] + DVi2[i]; 

//        fprintf(stderr,"N: %lf \n",N[i]);


          fprintf(fptr,"%d,%f,%d,%f,N\n",t,vv,i,N[i]);
          fprintf(fptr,"%d,%f,%d,%f,S\n",t,vv,i,S[i]);
          fprintf(fptr,"%d,%f,%d,%f,I1\n",t,vv,i,I1[i]);
          fprintf(fptr,"%d,%f,%d,%f,I2\n",t,vv,i,I2[i]);
          fprintf(fptr,"%d,%f,%d,%f,R1\n",t,vv,i,R1[i]);
          fprintf(fptr,"%d,%f,%d,%f,R2\n",t,vv,i,R2[i]);
          fprintf(fptr,"%d,%f,%d,%f,V\n",t,vv,i,V[i]);
          fprintf(fptr,"%d,%f,%d,%f,XIVI1\n",t,vv,i,XIV1[i]);
          fprintf(fptr,"%d,%f,%d,%f,XIV2\n",t,vv,i,XIV2[i]);
          fprintf(fptr,"%d,%f,%d,%f,Xsi1\n",t,vv,i,Xsi1[i]);
          fprintf(fptr,"%d,%f,%d,%f,Xsi2\n",t,vv,i,Xsi2[i]);
          fprintf(fptr,"%d,%f,%d,%f,XD\n",t,vv,i,XD[i]);
          fprintf(fptr,"%d,%f,%d,%f,VR1\n",t,vv,i,VR1[i]);
 //         fprintf(fptr,"%d,%f,%d,%f,VR2\n",t,vv,i,VR2[i]);
          fprintf(fptr,"%d,%f,%d,%f,H1\n",t,vv,i,H1[i]);
 //         fprintf(fptr,"%d,%f,%d,%f,H2\n",t,vv,i,H2[i]);

        }

        
 //        fprintf(stderr,"TOTAL MORT: %f, TOTAL BIRTHS: %f \n",total(D),total(SB));
 //        fflush(stderr);

        total_lambda = 0;

        for(int c = 0; c < AGES; c++){
            Xsi1[c] = 0;
            Xsi2[c] = 0;
            Xr1[c] = 0;
            Y1[c] = 0;
            Y2[c] = 0;
            YV1[c] = 0;
            YV2[c] = 0;
            omega1[c] = 0;
            omega2[c] = 0;
            omega3[c] = 0;
            omega4[c] = 0;
            theta1[c] = 0;
            theta2[c] = 0;
            theta3[c] = 0;
            theta4[c] = 0;
            sv[c] = 0;
            Xsh1[c] = 0;
            Xsh2[c] = 0;
            Xr1i2[c] = 0;
            Xr2i1[c] = 0;
            Xvvi1[c] = 0;
            Xvvi2[c] = 0;
            Yh1r1[c] = 0;
            Yh2r2[c] = 0;
            Xvr1vi2[c] = 0;
            Xvr2vi1[c] = 0;
            Di1[c] = 0;
            DVi1[c] = 0;
            DVi2[c] = 0;
            Di2[c] = 0;
            Dh1[c] = 0;
            Dh2[c] = 0;
            Vh1[c] = 0;
            Vh2[c] = 0;
            Vr1[c] = 0;
            Vs[c] = 0;
            vh1vr1[c] = 0;
            vh2vr2[c] = 0;
            p_vs[c] = 0;
            r1v[c] = 0;
            r2v[c] = 0;
            IS[c] = 0;
            II1[c] = 0;
            IVI1[c] = 0;
            IR1[c] = 0;
            IH1[c] = 0;
            IV[c] = 0;
            II2[c] = 0; 
            IVI2[c] = 0;
            IR2[c] = 0;
            SD[c] = 0;
            SB[c] = 0;
            I1D[c] = 0;
            I2D[c] = 0;
            VI1D[c] = 0;
            VI2D[c] = 0;
            R2D[c] = 0;
            R1D[c] = 0;
            VR2D[c] = 0;
            VR1D[c] = 0;
            VD[c] = 0;
            H2D[c] = 0;
            H1D[c] = 0;
            IH2[c] = 0;
            D[c] = 0;
            IVR1[c] = 0;
            IVR2[c] = 0;
        }
        t += 1;
    }
   // fprintf(stderr,"FINISHED STOCH");
    //fflush(stderr);
    free(m);
    free(mu_i1);
    free(mu_i2);
    free(VC);
    free(theta);
    free(ICU_raio); 
    free(psi);
    free(school_time);
    free(im_prop);
    free(S);
    free(I1);
    free(VI1);
    free(R1);
    free(VR1);
    free(VR2);
    free(H1);
    free(V);
    free(I2);
    free(VI2);
    free(R2);
    free(H2);
    free(Xsi1);
    free(Xsi2);
    free(XIV1);
    free(XIV2);
    free(XI1);
    free(XI2);
    free(XD);
    for(int i = 0;  i < AGES; i++){
        free(cm_overall[i]);
    //    free(cm_school[i]);
        free(M[i]);
    }
    free(cm_school);
    free(cm_overall);
    free(M);
    return; 
}
   //fprintf(stderr,"STARTING TO ENTER LOOP!");
   //fflush(stderr);
//    for(int c = 0; c < AGES; c++){
//        fprintf(stderr,"VI VALS: %f",VI1[c]);
//        fflush(stderr);
//    }

    // set all transition arrays to zero


    //oupu becomes massive due o many years 
    //print only wha's needed 
    // 1) write a sprintf when loopintg over age, if 40 compartments, looping over age writing 40 comparmens
    // 2) no interesintg int everyhintg, intsead of printintg 40 comparmens, wrie 10 - someimes don' need oupu for everyday
    // Incidience of cases, print everyday 
    // Prevalnce do yearly 
    // Deahs should be daily, hospia
    // Want: intcidence of variables ha change int shor erm
    //    -> cases, deahs, hospilizaions
    // print out daily vaccintaion rae for debuggintg 
