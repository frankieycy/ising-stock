import numpy as np
import matplotlib.pyplot as plt
from matplotlib import rc

rc('text',usetex=True)
plt.figure(figsize=(14,4))

d = np.loadtxt("out/data.csv",delimiter=",",skiprows=1)
t = d[:,0] # time
p = d[:,1] # price
m = d[:,2] # magnitisation
u = d[:,3] # return

axes = plt.gca()
axes.set_ylim([-6,6]) # return

plt.plot(t,p,color='black')
plt.xlabel('time')
plt.ylabel('return')
plt.savefig('plt/plot.png',bbox_inches='tight',dpi=200)
