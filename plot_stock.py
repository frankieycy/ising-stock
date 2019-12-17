import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt

# generate plots and statistics with model outputs

mpl.rc('text',usetex=True)

loc   = 'plt/'
data  = 'ising' # ising, real or random walk

STOCK = 1 # flag for stock plots (visual)
STAT  = 0 # flag for stat analysis (analytical)

# ------------------------------------- #

if data == 'ising':
	f = 'out/data.csv' # ising stock
if data == 'real':
	f = 'out/hsi.csv' # real stock
if data == 'random walk':
	f = 'out/rw.csv' # random walk

# 0th-col: time; 1st-col: price
d = np.loadtxt(f,delimiter=',',skiprows=1)
t = d[:-1,0] # time
p = d[:,1]   # price
u = (p[1:]-p[:-1])/p[:-1] # return

if STOCK == 1:
	# normalise return
	mu0 = np.mean(u)
	sig0 = np.std(u)
	u_norm = (u-mu0)/sig0

	# return: skewness & kurtosis
	from scipy.stats import skew,kurtosis
	print('skew:',round(skew(u_norm),4))
	print('kurt:',round(kurtosis(u_norm),4))

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

"""
	# return: 1-day lag plot
	fig = plt.figure(figsize=(6,6))
	ax = plt.gca()
	ax.set_xlim([-10,10])
	ax.set_ylim([-10,10])
	plt.scatter(u_norm[:-1],u_norm[1:],color='black',marker='.')
	plt.title('Lagged return: '+data)
	plt.xlabel('return today')
	plt.ylabel('return tomorrow')
	fig.tight_layout()
	fig.savefig('return.png')
"""

# ------------------------------------- #

def autocor(ts,tau):
	# autocorrelation for time series with lag tau
	T = len(ts)
	mu = np.mean(ts)
	return sum((ts[0:T-tau]-mu)*(ts[tau:T]-mu))/sum((ts-mu)**2)

if STAT == 1:
	taus = range(1,500) # time lags
	rho = [] # autocor of return
	rho_a = [] # autocor of abs return
	for tau in taus:
		rho.append(autocor(u,tau))
		rho_a.append(autocor(np.abs(u),tau))

	# power-law fit; slope = alpha
	log_taus = []
	log_rho_a = []
	for i in range(len(rho_a)):
		if rho_a[i]>0:
			log_taus.append(np.log(taus[i]))
			log_rho_a.append(np.log(rho_a[i]))

	alpha = -np.polyfit(log_taus,log_rho_a,deg=1)[0]
	print('alpha:',round(alpha,4))
"""
	# return: autocorrelation
	fig = plt.figure(figsize=(6,4))
	plt.plot(rho,color='black')
	plt.axhline(y=+1.96/np.sqrt(len(u)),color='black',linestyle='--')
	plt.axhline(y=-1.96/np.sqrt(len(u)),color='black',linestyle='--')
	plt.title('Autocorrelation of return: '+data)
	plt.xlabel('lag')
	plt.grid()
	fig.tight_layout()
	fig.savefig(loc+data+'_autocor.png')

	# abs return: autocorrelation
	fig = plt.figure(figsize=(6,4))
	plt.plot(rho_a,color='black')
	plt.title('Autocorrelation of absolute return: '+data)
	plt.xlabel('lag')
	plt.grid()
	fig.tight_layout()
	fig.savefig(loc+data+'_autocor-abs.png')
"""
