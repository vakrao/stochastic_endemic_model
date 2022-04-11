import os
import sys
import glob
import string
import pandas as pd

#finds total sum of all variables regardless of age
vax_levels = 4 
starting_config =0 
total_mean_name = "all_vax_mean_values_exp1.csv"
total_quantile_name = "all_vax_quantile_values_exp1.csv"
titles = []
qtitles = []
#finds mean and quantile for each vaccine level 
for i in range(starting_config,vax_levels):
    file_name = "total_sim_vax_level_"+str(i)+"0.csv"
    raw_data = pd.read_csv(file_name)
    mean_aged_data = raw_data.groupby(['t','vartype','age']).mean()
    aged_quantile_data = raw_data.groupby(['t','vartype','age']).quantile(q=[0.25,0.5,0.75])
    raw_data = raw_data.drop("age",1)
    mean_data = raw_data.groupby(['t','vartype']).mean()
    quantile_data = raw_data
    quantile_data["quantile"] = raw_data.groupby(['t','vartype']).quantile(q=[0.25,0.5,0.75])

    mean_name = "total_sim_vax_level_mean_"+str(i)+"0.csv"
    mean_age_name = "total_sim_vax_ages_level_mean_"+str(i)+"0.csv"
    quantile_name = "total_sim_vax_level_quantile_"+str(i)+"0.csv"
    quantile_age_name = "total_sim_vax_ages_level_quantile_"+str(i)+"0.csv"
    mean_aged_data.to_csv(mean_age_name)
    mean_data.to_csv(mean_name)
    quantile_data.to_csv(quantile_name)
    aged_quantile_data.to_csv(quantile_age_name)
    print("FINISHED FOR: %d",i)

#joins all values together
for k in range(starting_config,vax_levels):
    string_val = str(k)
    title_name = "total_sim_vax_level_mean_"+string_val+"0.csv"
    quantile_title_name = "total_sim_vax_level_quantile_"+string_val+"0.csv"
    titles.append(title_name)
    qtitles.append(quantile_title_name)

combined_mean_csv = pd.concat([pd.read_csv(f) for f in titles])
combined_mean_csv.to_csv( total_mean_name, index=False, encoding='utf-8-sig')
combined_quantile_csv = pd.concat([pd.read_csv(f) for f in qtitles])
combined_quantile_csv.to_csv( total_quantile_name, index=False, encoding='utf-8-sig')
