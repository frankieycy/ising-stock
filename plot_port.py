import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.figure import Figure

mpl.rc('text',usetex=True)
Figure(tight_layout=True)

loc = 'plt/'
mod = 'poisson' # base or poisson

f = 'out/port.csv'
d = np.loadtxt(f,delimiter=',',skiprows=1)
v = d[:,1]

# portfolio: distribution
fig = plt.figure(figsize=(6,4))
plt.hist(v,color='black',bins=25,density=True)
plt.title('Distribution of portfolio value: '+mod)
plt.xlabel('portfolio value')
fig.tight_layout()
fig.savefig(loc+mod+'_port.png')
