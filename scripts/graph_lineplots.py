import pandas as pd 
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt

sim_data = pd.read_csv("total_mean_values.csv")
#wide_sim_data = pd.read_csv("/results/sim_one.csv")
wide_sim_quantile_data = pd.read_csv("total_sim_vax_level_quantile_10.csv")
infec_title = "LINC_exp1.png"
death_title = "LDTH_exp1.png"
#print("I1 : ",max(sim_data["I1"]))
#print("I2 : ",max(sim_data["I2"]))
#print("XSi1 : ",max(sim_data["XSi1"]))
#print("XSi2: ",max(sim_data["XSi2"]))
#print("VI1 : ",max(sim_data["VI1"]))
#print("VI2 : ",max(sim_data["VI2"]))
#print("N : ",max(sim_data["N"]))
valueX1 = "Xsi1"
valueX2 = "Xsi2"
valueD = "XD"
print(sim_data['vartype'])
infected_data = sim_data.query('vartype == @valueX1')
infected_data["rolling"] =infected_data["value"].rolling(window=7).sum()
death_data = sim_data.query('vartype == @valueD')
death_data["rolling"] =death_data["value"].rolling(window=7).sum()
#sns.lineplot(data=infected_data,x="t",y="value",hue="vartype")
g = sns.FacetGrid(infected_data,col="vv",height=2.5,col_wrap=3)
g2 = sns.FacetGrid(death_data,col="vv",height=2.5,col_wrap=3)

g.map(sns.lineplot,"t","rolling")
g.set(ylim=(0,20000))
g2.map(sns.lineplot,"t","rolling")
#g2.set(ylim=(0,20000))
fig = g.savefig(infec_title)
fig2 = g2.savefig(death_title)
plt.show()
