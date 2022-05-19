#!/bin/bash

#SBATCH -J stoch_model_run #job name                                      
#SBATCH -p general #submit to the general partition
#SBATCH -o filename_%j.txt #output file name (%j is replaced by job ID)
#SBATCH -e filename_%j.err #standard error file name (%j is replaced by job ID)
#SBATCH --mail-type=ALL #send email for ALL, BEGIN, END, or FAIL
#SBATCH --mail-user=vakrao@iu.edu #the email address for job-related emails
#SBATCH --nodes=1 #number of requested nodes
#SBATCH --cpus-per-task=1
#SBATCH --mem=100gb #memory per node
#SBATCH --ntasks-per-node=1 #number of tasks per node
#SBATCH --time=10:00:00 #requested time


#load modules
module load gsl
#module load openmp
module swap gcc/6.3.0 gcc/10.2.0
modue load python

#path to your working directory
#cd full/path/to/your/directory

#commands to execute in following lines
echo '#### Move to source directory ####'
#cd /N/slate/vakrao/stoch_endemic_model/src
echo $PWD
echo '#### Running model! ####' 
#gcc -fopenmp -O2 -mtune=native -march=native -Wall -c *.c 
#gcc -o model *.o -lgsl -lgslcblas -lm -fopenmp
gcc -c *.c 
gcc -o model *.o -lgsl -lgslcblas -lm 
#export OMP_NUM_THREADS=12
srun model 10 0 
srun model 10 10 
srun model 10 20 
srun model 10 30 
srun model 10 40 
srun model 10 50 
srun model 10 60 
srun model 10 70 
srun model 10 80 
srun model 10 90 
srun model 10 100 
echo '### Moving files to another folder !'
mv *.csv ../scripts/
cd ../scripts/
module load python
#echo '###### Starting paste!######'
python3 seq_csv.py
#echo '###### Finished paste!######'
#echo '###### Starting mean!######'
python3 find_mean.py
#echo '###### Finished mean!######'
zip -9 exp1_quantile total_sim_vax_level_quantile_*.csv
zip -9 exp1_mean total_sim_vax_level_mean_*.csv
#zip -9 exp1_ages_mean total_sim_vax_ages_level_mean_*.csv
#zip -9 exp1_ages_quantile total_sim_vax_ages_level_quantile_*.csv
echo '#### Finished zipping! ####'
python graph_lineplots.py
#echo '#### Finished plotting! ####'
#
