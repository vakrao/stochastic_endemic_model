#!/bin/sh

#SBATCH -J Viz_run #job name                                      
#SBATCH -p general #submit to the general partition
#SBATCH --mail-type=ALL #send email for ALL, BEGIN, END, or FAIL
#SBATCH --mail-user=shedmund@iu.edu #the email address for job-related emails
#SBATCH --nodes=1 #number of requested nodes
#SBATCH --cpus-per-task=8
#SBATCH --mem=100gb #memory per node
#SBATCH --ntasks-per-node=1 #number of tasks per node
#SBATCH --time=4:00:00 #requested time

module load r

R CMD BATCH viz_script.R
#R CMD BATCH R0_3.R
R CMD BATCH R0_5.R
R CMD BATCH R0_7.R
R CMD BATCH R0_9.R

#R CMD BATCH 100_days.R
#R CMD BATCH 150_days.R

#R CMD BATCH VE_30.R
#R CMD BATCH VE_40.R
