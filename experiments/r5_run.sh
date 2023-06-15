#!/bin/bash

#SBATCH -J 50_r5 #job name                                      
#SBATCH -p general #submit to the general partition
#SBATCH --mail-type=ALL #send email for ALL, BEGIN, END, or FAIL
#SBATCH --mail-user=shedmund@iu.edu #the email address for job-related emails
#SBATCH --error=%j.err
#SBATCH --nodes=1 #number of requested nodes
#SBATCH --cpus-per-task=11
#SBATCH --mem=150gb #memory per node
#SBATCH --ntasks-per-node=1 #number of tasks per node
#SBATCH --time=12:00:00 #requested time

# R0, gamma, duration
export RUNS=50
export FOLDER="R0_5"
export FOLDERDIR="R0"
export BIGFOLDER="/N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw/"
export INPUT_PARAMS="/N/project/endemic_covid/experiments/specs/$FOLDER.csv"
export OUTPUT_PARAMS="/N/project/endemic_covid/data/specs/$FOLDER.csv"
export TYPE="age"

#load modules
module load gsl
module swap gcc/6.3.0 gcc/10.2.0
module load python

#path to your working directory
#cd full/path/to/your/directory

#Now, we compile all C files and create the object file...
#gcc -c *.c 
echo '#### Compiling model! ####' 
#gcc -o model *.o -lgsl -lgslcblas -lm 
echo '#### Running model! ####' 
#Running model for different vaccine configs, 50 runs
srun model $RUNS 00 $AGE $PARAMS $BIGFOLDER  
srun model $RUNS 10 $AGE $PARAMS $BIGFOLDER 
srun model $RUNS 20 $AGE $PARAMS $BIGFOLDER 
srun model $RUNS 30 $AGE $PARAMS $BIGFOLDER 
srun model $RUNS 40 $AGE $PARAMS $BIGFOLDER  
srun model $RUNS 50 $AGE $PARAMS $BIGFOLDER 
srun model $RUNS 60 $AGE $PARAMS $BIGFOLDER 
srun model $RUNS 70 $AGE $PARAMS $BIGFOLDER 
srun model $RUNS 80 $AGE $PARAMS $BIGFOLDER 
srun model $RUNS 90 $AGE $PARAMS $BIGFOLDER  
srun model $RUNS 100 $AGE $PARAMS $BIGFOLDER   
rm *.o
#Now, we move all files to the temp folder!
echo '### Moving files to another folder !'
cp $INPUT_PARAMS $OUTPUT_PARAMS

cd /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw

for i in {0..10}; do cat age_run_*_"$i"0.csv > /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/total/total_"$i"0.csv; done
#mv age_run_*_*.csv /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw

module load r

R CMD BATCH /N/project/endemic_covid/scripts/viz_script.R
R CMD BATCH /N/project/endemic_covid/scripts/R0_5.R
