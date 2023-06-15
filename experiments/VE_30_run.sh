#!/bin/bash

#SBATCH -J 50_ve_30 #job name                                      
#SBATCH -p general #submit to the general partition
#SBATCH --mail-type=ALL #send email for ALL, BEGIN, END, or FAIL
#SBATCH --mail-user=shedmund@iu.edu #the email address for job-related emails
#SBATCH --error=%j.err
#SBATCH --nodes=1 #number of requested nodes
#SBATCH --cpus-per-task=11
#SBATCH --mem=100gb #memory per node
#SBATCH --ntasks-per-node=1 #number of tasks per node
#SBATCH --time=8:00:00 #requested time

export RUNS=50
export FOLDER="30"
export FOLDERDIR="VE"
export BIGFOLDER="/N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw/"
export PARAMS="VE_30.csv"

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

cd /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw

for i in {0..10}; do cat age_run_*_"$i"0.csv > /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/total/total_"$i"0.csv; done
#mv age_run_*_*.csv /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw

cd /N/project/endemic_covid/scripts/

module load r

R CMD BATCH viz_script.R
R CMD BATCH VE_30.R
