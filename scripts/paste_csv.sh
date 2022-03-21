#!/bin/bash
#SBATCH --mail-user=vakrao@iu.edu
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1
#SBATCH --time=0-10:00:00
#SBATCH --partition=general
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --job-name=csv_paste_job
#SBATCH --mem=200gb #memory per node
#SBATCH --output=log

######  Module commands #####
module load python/3.6.9



######  Job commands go below this line #####
echo '###### move to script dir ######'
cd ~/stoch_endemic_model/scripts
echo '###### Directory Changed! ######'

echo '###### Running csv script ######'
python3 join_csv.py
echo '###### Run Complete! ######'
