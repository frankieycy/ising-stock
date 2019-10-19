import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

mpl.rc('text',usetex=True)

loc = 'plt/'
data = 'ising' # ising or real

if data == 'ising':
	f = 'out/data.csv' # ising stock
if data == 'real':
	f = '0001.HK.csv' # real stock

d = np.loadtxt(f,delimiter=',',skiprows=1)
t = d[:-1,0] # time
p = d[:,1]   # price
u = (p[1:]-p[:-1])/p[:-1] # return

# normalise return
mu0 = np.mean(u)
sig0 = np.std(u)
u_norm = (u-mu0)/sig0

# return: time series
fig = plt.figure(figsize=(10,4))
ax = plt.gca()
ax.set_ylim([-6,6])
plt.plot(t,u_norm,color='black')
plt.title('Stock return: '+data)
plt.xlabel('time')
plt.ylabel('return')
fig.tight_layout()
fig.savefig(loc+data+'_return.png')

# price: time series
fig = plt.figure(figsize=(10,4))
plt.plot(t,p[:-1],color='black')
plt.title('Stock price: '+data)
plt.xlabel('time')
plt.ylabel('price')
fig.tight_layout()
fig.savefig(loc+data+'_price.png')

# return: quantile plot
n = np.random.normal(size=len(t))
x = np.array([np.min(n),np.max(n)])
y = x
fig = plt.figure(figsize=(6,4))
plt.scatter(np.sort(n),np.sort(u_norm),color='black',marker='.')
plt.plot(x,y,color='black') # theoretical fit
plt.title('Quantile plot: '+data)
plt.xlabel('theoretical normal')
plt.ylabel('sample')
fig.tight_layout()
fig.savefig(loc+data+'_qq.png')

# return: distribution
fig = plt.figure(figsize=(6,4))
_,bins,_ = plt.hist(u_norm,color='black',bins=25,density=True,range=[-4,4])
plt.plot(bins,1/np.sqrt(2*np.pi)*np.exp(-bins**2/2),color='grey')
plt.title('Distribution of return: '+data)
plt.xlabel('return')
fig.tight_layout()
fig.savefig(loc+data+'_dist.png')
