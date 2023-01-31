#!/bin/bash

#SBATCH -J cat_files #job name                                      
#SBATCH -p general #submit to the general partition
#SBATCH --mail-type=ALL #send email for ALL, BEGIN, END, or FAIL
#SBATCH --mail-user=shedmund@iu.edu #the email address for job-related emails
#SBATCH --error=%j.err
#SBATCH --nodes=1 #number of requested nodes
#SBATCH --cpus-per-task=4
#SBATCH --mem=100gb #memory per node
#SBATCH --ntasks-per-node=1 #number of tasks per node
#SBATCH --time=1:00:00 #requested time

export FOLDER="40"
export FOLDERDIR="VE"
#load modules
module load gsl
module swap gcc/6.3.0 gcc/10.2.0

cd /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw
echo '###### Creating large files! ######'
for i in {0..10}; do cat age_run_*_"$i"0.csv > /N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/total/total_"$i"0.csv; done
