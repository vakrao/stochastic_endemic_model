import os
import sys
import glob
import string
import pandas as pd
#os.chdir("../src")
extension = 'csv'
titles = []
total_mean_name = "total_mean_vax_levels.csv"
for k in range(0,11):
    string_val = str(k)
    title_name = "../results/exp_1/raw_runs/run_*_"+string_val+"0.csv"
    print(title_name)
    all_filenames = [i for i in glob.glob(title_name)]
    #combine all files in the list
    for j in range(0,11):
        all_filenames[j] = "../results/exp_1/raw_runs/"+os.path.basename(all_filenames[j])
    print(all_filenames)
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ])
    #export to csv
    total_sim_string_name = "total_sim_vax_level_"+str(k)+"0.csv"
    print("finished one")
    combined_csv.to_csv( total_sim_string_name, index=False, encoding='utf-8-sig')

for k in range(0,11):
    string_val = str(k)
    title_name = "total_sim_vax_level_mean_"+string_val+"0.csv"
    titles.append(title_name)

combined_csv = pd.concat([pd.read_csv(f) for f in titles])
combined_csv.to_csv( total_mean_name, index=False, encoding='utf-8-sig')