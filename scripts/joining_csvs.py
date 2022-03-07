import os
import sys
import glob
import string
import pandas as pd
#os.chdir("../src")
extension = 'csv'
for k in range(0,10):
    string_val = str(k)
    title_name = "../src/run_*_"+string_val+"0.csv"
    all_filenames = [i for i in glob.glob(title_name)]
    #combine all files in the list
    for j in range(0,11):
        all_filenames[j] = "../src/"+os.path.basename(all_filenames[j])
    combined_csv = pd.concat([pd.read_csv(f) for f in all_filenames ])
    #export to csv
    total_sim_string_name = "total_sim_vax_level_"+str(k)+"0.csv"
    combined_csv.to_csv( total_sim_string_name, index=False, encoding='utf-8-sig')