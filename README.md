## Ising Stock Market

* PHYS4061 Project B: Monte Carlo Simulation
* Last: 18/10/2019

## Model

* market: (i) fundamentalists (implicit), (ii) interacting traders
* interacting traders placed in 2d lattice with periodic boundary (2d Ising)
* spin (+/-): buy/sell decision
* _locally_ seeking ferromagnetic order (nbor-interaction)
* indiv interacting traders influenced by nbor
* _globally_ escaping ferromagnetic order (external field)
* returns & stock prices updated by net magnetisation
* snaphots taken per 1000 time steps (black & white as opposite spins):

	![](ani.gif)

## Params

## Codes

* folders:
	- `out/`: data outputs
	- `plt/`: python plots
	- `bak/`: backups
	- `rw/`: random walk model (later)
* outputs:
	- data.csv: price, magnetisation, return as time series (ts)
	- spin.csv: spins on lattice
	- port.csv: portfolio values of traders

## Results

* plots:
	- stock price as ts
	- return as ts: volatility clustering
	- distribution of return: non-normal
	- quantile plot of return: heavy tails
	- autocorrelation (later)

## References

* Ising Model in Finance (Pavel Dvorak)
* Dynamics of price and trading volume in a spin model of stock markets with heterogeneous agents (Bornholdt)
