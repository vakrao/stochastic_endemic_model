import pandas as pd 
import seaborn as sns 
import numpy as np 
import matplotlib.pyplot as plt

sim_data = pd.read_csv("../results/total_all_comps_sim.csv")
wide_sim_data = pd.read_csv("../results/sim_one.csv")
#print("I1 : ",max(sim_data["I1"]))
#print("I2 : ",max(sim_data["I2"]))
#print("XSi1 : ",max(sim_data["XSi1"]))
#print("XSi2: ",max(sim_data["XSi2"]))
#print("VI1 : ",max(sim_data["VI1"]))
#print("VI2 : ",max(sim_data["VI2"]))
#print("N : ",max(sim_data["N"]))
#infected_data = wide_sim_data.query('vartype == "VI1 \n1"')
sns.lineplot(data=sim_data,x="t",y="value",hue="vartype")
plt.show()
