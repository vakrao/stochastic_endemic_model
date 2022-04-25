import pandas as pd 
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt

infec_limit = 20000
death_limit = 20000

#sim_data = pd.read_csv("all_vax_mean_values_exp1.csv")
sim_data = pd.read_csv("run_0_00.csv")
sim_data["vv"] = sim_data["vv"].round(2)
sim_data["value"] = sim_data["value"].round(2)
#wide_sim_data = pd.read_csv("/results/sim_one.csv")
#wide_sim_quantile_data = pd.read_csv("total_mean_vax_levels.csv")
valueX1 = "Xsi1"
valueX2 = "Xsi2"
valueD = "XD"
infected_data = sim_data.query('vartype == @valueX1')
death_data = sim_data.query('vartype == @valueD')
sim_data["vv"].round(2)
#sns.lineplot(data=infected_data,x="t",y="value",hue="vartype")
g = sns.FacetGrid(infected_data,col="vv",height=2.5,col_wrap=3)
g.map(sns.lineplot,"t","value")
g2 = sns.FacetGrid(death_data,col="vv",height=2.5,col_wrap=3)
g2.map(sns.lineplot,"t","value")
#g.set(ylim=(0,infec_limit))
#g2.set(ylim=(0,death_limit))
fig = g.savefig("exp1_meaninfec_100.png")
fig2 = g2.savefig("exp1_meandeath_100.png")

