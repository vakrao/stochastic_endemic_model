#!/bin/bash

#SBATCH -J 20_endemic_run #job name                                      
#SBATCH -p general #submit to the general partition
#SBATCH --mail-type=ALL #send email for ALL, BEGIN, END, or FAIL
#SBATCH --mail-user=vakrao@iu.edu #the email address for job-related emails
#SBATCH --error=%j.err
#SBATCH --nodes=1 #number of requested nodes
#SBATCH --cpus-per-task=5
#SBATCH --mem=100gb #memory per node
#SBATCH --ntasks-per-node=1 #number of tasks per node
#SBATCH --time=5:00:00 #requested time

export PARAMS="one_strain.csv"
export RUNS=20
export FOLDER="R0_5"
export FOLDERDIR="R0"
export BIGFOLDER="/N/project/endemic_covid/data/raw/R0/R0_5/raw/"
#load modules
module load gsl
module swap gcc/6.3.0 gcc/10.2.0
module load python

#path to your working directory
#cd full/path/to/your/directory

#Now, we compile all C files and create the object file...
gcc -c *.c 
echo '#### Compiling model! ####' 
gcc -o model *.o -lgsl -lgslcblas -lm 
echo '#### Running model! ####' 
#Running model for different vaccine configs, 50 runs
srun model $RUNS 00 0 $PARAMS $BIGFOLDER 
srun model $RUNS 10 0 $PARAMS $BIGFOLDER 
srun model $RUNS 20 0 $PARAMS $BIGFOLDER 
srun model $RUNS 30 0 $PARAMS $BIGFOLDER 
srun model $RUNS 40 0 $PARAMS $BIGFOLDER  
srun model $RUNS 50 0 $PARAMS $BIGFOLDER 
srun model $RUNS 60 0 $PARAMS $BIGFOLDER 
srun model $RUNS 70 0 $PARAMS $BIGFOLDER 
srun model $RUNS 80 0 $PARAMS $BIGFOLDER 
srun model $RUNS 90 0 $PARAMS $BIGFOLDER  
srun model $RUNS 100 0 $PARAMS $BIGFOLDER   
rm *.o
#Now, we move all files to the temp folder!
echo '### Moving files to another folder !'
cp one_strain.csv /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER
#mv age_run_*_*.csv /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw

cd /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw
