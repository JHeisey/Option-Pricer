import random
import numpy as np
import csv
np.set_printoptions(suppress=True)
np.set_printoptions(precision=4)

rows = 1000000

Smin = 90
Smax = 91
Kmin = 89
Kmax = 92
rmin = 0.05
rmax = 0.06
vmin = 0.20
vmax = 0.25

T = np.ones((rows,1))/52.0
S = np.random.uniform(low=Smin, high=Smax, size=(rows,1))
K = np.random.uniform(low=Kmin, high=Kmax, size=(rows,1))
r = np.random.uniform(low=rmin, high=rmax, size=(rows,1))
v = np.random.uniform(low=vmin, high=vmax, size=(rows,1))

vars = np.column_stack((S,K,r,v,T))

print(vars)

np.savetxt("optionPricing/variables.csv",vars,fmt='%1.4f', delimiter=",")
