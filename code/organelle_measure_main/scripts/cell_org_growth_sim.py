# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 15:03:42 2024

@author: Anang
"""

import numpy as np
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import poisson
import seaborn as sns


#%% Paramenters
'''
growth=y0+Ae^(x/t)

Ae^(growth*x)


poisson dist= ((lambda_value2**k)*np.exp(-lambda_value2))/math.factorial(k)


at measurement time their must be a poisson distribution of size of the cell

with this poisson distribution we must have a cell size,

this cell size function is nothing but the fuction of growth rate and 


'''



#get the poisson distribution for the size of the cell

lambda_param = 1

# Generate a random sample from a Poisson distribution
data = np.random.poisson(lambda_param, size=10000)
# Plot the histogram of the generated data
# plt.hist(data, bins=20, density=True, alpha=0.7, color='blue', edgecolor='black')
sns.kdeplot(data, bw=0.5)


plt.figure(1)
x_values = np.arange(0, 20)

x = np.arange(0, 30, 0.05)
y = poisson.pmf(x, mu=20, loc=0)*20
plt.scatter(x, y)

plt.figure(2)

sns.kdeplot(y, bw=0.5)

# showing the graph
plt.show()

#get the cell size function 
plt.figure(2)
v0=1
vcell1=v0*(pow(np.e, y/0.5))
vcell2=v0*(pow(np.e, y/0.4))
vcell3=v0*(pow(np.e, y/0.29))
vcell4=v0*(pow(np.e, y/0.16))
vcell5=v0*(pow(np.e, y/0.13))

plt.scatter(x, vcell1,label="vcell")
plt.scatter(x, vcell2,label="vcell")
plt.scatter(x, vcell3,label="vcell")
plt.scatter(x, vcell4,label="vcell")
plt.scatter(x, vcell5,label="vcell")

# showing the graph
plt.show()


fig, ax = plt.subplots()
ax.hist(vcell5)


#%%

f1=np.e(k)

kd=0.4
kf=0.3
kfu=0.1
ky=1.0-(kd+kf+kfu)

per=kd+kf
vac=kf+kfu
er=
golgi
mito=kf+kf
nonorg=

function=per+vac+er+golgi+mito+lD+nonorg
#%%

rc=1

ro=0.8

fraction=pow(ro,3)/pow(rc,3)





org0=0.6

org1=0.1

1/((0.6)**(1/3))

(0.6)**(1/3)

(0.84)**3