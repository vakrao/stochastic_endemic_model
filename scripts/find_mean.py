import os
import sys
import glob
import string
import pandas as pd

#finds total sum of all variables regardless of age
vax_levels = 10
for i in range(0,vax_levels):
    file_name = "total_sim_vax_level_"+str(i)+"0.csv"
    raw_data = pd.read_csv(file_name)
    mean_data = raw_data.groupby(['t','vartype']).mean()
    quantile_data = raw_data.groupby(['t','vartype']).quantile(q=[0.25,0.9])
    mean_name = "total_sim_vax_level_mean_"+str(i)+"0.csv"
    quantile_name = "total_sim_vax_level_quantile_"+str(i)+"0.csv"
    mean_data.to_csv(mean_name)
    quantile_data.to_csv(quantile_name)
