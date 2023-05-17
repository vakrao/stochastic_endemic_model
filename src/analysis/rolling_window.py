import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

dat = pd.read_csv('total_100.csv')
wsize = 365
total_inc_data = dat[dat.vartype == "Xsi1"]
total_inc_data['value'] = total_inc_data['value'].astype(float)
total_inc_data['t'] = total_inc_data['t'].astype(int)
print(total_inc_data['t'])

total_inc_data['year'] = (total_inc_data["t"])/365
total_inc_data = total_inc_data.groupby(["t","vartype","year"]).sum()
total_inc_data = total_inc_data.rolling(wsize).mean(wsize)
total_inc_data = total_inc_data.apply(np.log10)

sns.lineplot(data=total_inc_data,x='year',y="value")
