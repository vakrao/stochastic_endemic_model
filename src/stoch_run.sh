#!/bin/bash

#SBATCH -J 50_endemic_run #job name                                      
#SBATCH -p general #submit to the general partition
#SBATCH -o endemic_%j.txt #output file name (%j is replaced by job ID)
#SBATCH -e endemic_%j.err #standard error file name (%j is replaced by job ID)
#SBATCH --mail-type=ALL #send email for ALL, BEGIN, END, or FAIL
#SBATCH --mail-user=vakrao@iu.edu #the email address for job-related emails
#SBATCH --nodes=1 #number of requested nodes
#SBATCH --cpus-per-task=5
#SBATCH --mem=100gb #memory per node
#SBATCH --ntasks-per-node=1 #number of tasks per node
#SBATCH --time=10:00:00 #requested time

export PARAMS = "one_strain.csv"
export RUNS = 50
export FOLDER = "50_run"

#load modules
module load gsl
module swap gcc/6.3.0 gcc/10.2.0
modue load python

#path to your working directory
#cd full/path/to/your/directory

#Now, we compile all C files and create the object file...
gcc -c *.c 
echo '#### Compiling model! ####' 
gcc -o model *.o -lgsl -lgslcblas -lm 
echo '#### Running model! ####' 
#Running model for different vaccine configs, 50 runs
srun model $RUNS 0 0 $PARAMS  
srun model $RUNS 10 0 $PARAMS  
srun model $RUNS 20 0 $PARAMS  
srun model $RUNS 30 0 $PARAMS  
srun model $RUNS 40 0 $PARAMS   
srun model $RUNS 50 0 $PARAMS  
srun model $RUNS 60 0 $PARAMS  
srun model $RUNS 70 0 $PARAMS  
srun model $RUNS 80 0 $PARAMS  
srun model $RUNS 90 0 $PARAMS   
srun model $RUNS 100 0 $PARAMS   
rm *.o
#Now, we move all files to the scripts folder!
echo '### Moving files to another folder !'
mv *.csv ../scripts/
cd ../scripts/
module load python
echo '###### Creating large files! ######'
python3 seq_csv.py
mv total_vax_*.csv /N/slate/project/endemic_covid/data/raw/$FOLDER

#echo '###### Finished paste!######'
#echo '###### Starting mean!######'
#python3 age_mean.py
#echo '###### Finished mean!######'
#zip -9 exp1_quantile total_sim_vax_level_quantile_*.csv

#zip -9 exp1_mean total_sim_vax_level_mean_*.csv
#zip -9 exp1_ages_mean total_sim_vax_ages_level_mean_*.csv
#zip -9 exp1_ages_quantile total_sim_vax_ages_level_quantile_*.csv
#echo '#### Finished zipping! ####'
#python graph_lineplots.py
echo '#### Finished plotting! ####'

