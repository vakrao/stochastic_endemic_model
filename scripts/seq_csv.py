import os
import sys
import glob
import string
import pandas as pd
import multiprocessing 

#os.chdir("../src")
#extension = 'csv'

if __name__=="__main__":
    print("hello!")
    total_mean_name = "total_mean_vax_levels.csv"
    total_quantile_name = "total_quantile_vax_levels.csv"
    n = 1
    starting_config = 0
    for k in range(starting_config,n):
        string_val = str(k)
        title_name = "run_*_"+string_val+"0.csv"
        all_filenames = [i for i in glob.glob(title_name)]
        title_name = "run_*_"+string_val+"0.csv"
        all_filenames = [i for i in glob.glob(title_name)]
        #combine all files in the list
        for j in range(0,len(all_filenames)):
            all_filenames[j] = os.path.basename(all_filenames[j])
        combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ])
        print("combined!!")
        name = "total_sim_vax_level_"+str(k)+"0.csv"
        combined_csv.to_csv(name,index=False,encoding='utf-8-sig')
        
