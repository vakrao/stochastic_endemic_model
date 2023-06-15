import numpy as np
import pandas as pd
import csv

fname = "../params/overall_contacts.csv"
#fname = "../params/school_contacts.csv"
full_conts = np.loadtxt(fname,delimiter=",",dtype=float)

new_bins = 17
old_bins = 85
increment = 5
new_M = []
M = []

for i in range(0,new_bins):
    new_M.insert(i,[])
    for j in range(0,new_bins):
        new_M[i].append(0)

for i in range(0,new_bins):
    sr = (i)*increment
    er = sr + increment
    for j in range(0,new_bins):
        sc = (j)*increment
        ec = sc + increment
        new_M[i][j] = np.sum(full_conts[sr:er,sc:ec])

#for i in range(0,old_bins):
#    M.insert(i,[])
#    for j in range(0,old_bins):
#        M[i].append(0)
#
#for i in range(0,old_bins):
#    for j in range(0,old_bins):
#        r = int(i/new_bins)
#        c = int(j/new_bins)
#        M[i][j] = new_M[r][c]
#
        
if fname == "../params/overall_contacts.csv":
    with open('../params/overall_17_contacts.csv','w',newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(new_M)
if fname == "../params/school_contacts.csv":
    with open('../params/school_17_contacts.csv','w',newline='') as csvfile:
        writer = csv.writer(csvfile)
        writer.writerows(new_M)






