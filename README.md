# Stochastic Endemic Model
This repository provides the source code for a realistic age-structured epidemic model, where we compute the incidence of SARS-CoV2 cases over a span of 
20 years. Using this model, we are interested in testing the effect that different vaccine coverage schemes will have on disease incidence. 

## Getting started
### Prerequistes
GSL version 2.6 must be installed 
GCC 10.2.0
R 4.0+ (Tidyverse) 
		
## Running Model
To run the model, you first must compile and generate the executable file. 
We have provided a shell script to compile the executable 
After running the executable, move the executable to the experiments folder
Then, either create a custom SLURM job script or run one of the existing scripts 
`cd /src
 sh run.sh 
cd ../experiments/
sbatch r5_run.script ' 

## Creating a custom SLURM script 

The model currently has 5 input parameters 
	
1. Number of runs. Larger amount will generate more stochastic simulations.
2. Vaccine coverage level. (00,10,20,30,40,50,60,70,80,90,100)
3. Model type. "total, "age", or "category"
- Total saves totals of all model states
- Age saves individuasl values of each age group for each time step (massive!) 
- Category categorizes ages into 6 differnet age groups (0-4,5-12,13-17,18-49,50-64,65+) and saves N, Xsi, and XV 
4. Parameter .csv file location
5. Folder specifying where to save files. 
	
`srun $RUNS $COVERAGE $MODELTYPE $PARAMETERS $FOLDER`

## Acknowledgements
- Co-collaborator: Seth Edmunds
- PIs: Marco Ajelli, Maria Litvinova 
- Thank you to the MOBS lab for providing open-source contact network data
- Mortality, demographic data publicly available from CDC 
	  
		
