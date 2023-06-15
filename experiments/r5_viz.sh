#!/bin/bash
#SBATCH -J r5_viz #job name                                      
#SBATCH -p general #submit to the general partition
#SBATCH --mail-type=ALL #send email for ALL, BEGIN, END, or FAIL
#SBATCH --mail-user=vakrao@iu.edu #the email address for job-related emails
#SBATCH --error=%j.err
#SBATCH --nodes=1 #number of requested nodes
#SBATCH --cpus-per-task=11
#SBATCH --mem=150gb #memory per node
#SBATCH --ntasks-per-node=1 #number of tasks per node
#SBATCH --time=12:00:00 #requested time

export RUNS=50
export FOLDER="R0_5"
export FOLDERDIR="R0"
export BIGFOLDER="/N/project/endemic_covid/data/raw/$FOLDERDIR/$FOLDER/raw/"
export PARAMS="$FOLDER.csv"

#load modules
module load gsl
module swap gcc/6.3.0 gcc/10.2.0
module load python

#path to your working directory
#cd full/path/to/your/directory
R CMD BATCH /N/project/endemic_covid/scripts/viz_script.R
R CMD BATCH /N/project/endemic_covid/scripts/R0_5.R
