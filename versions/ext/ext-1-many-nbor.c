#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tool.c"
#include "vector.c"

/*
# Model
many-nbor but weak-coupling model
traders subject to more diverse but weaker local couplings

# Definitions
T     total time
dt    num of iter in each MC step (each trade day)
t     current time (day)

N     num of traders
K     cutoff distance (to pick nbor)
R     lattice coord (vec[N][N])
nb    num of nbor
Nb    nbor list (int[N][N])
D     pbc nbor distances (float[N][N])
V     value of one's portfolio (initial=0)

A     mean field coupling strength
B     transition probability control
J     coupling strength

S     spin (+1 or -1) (float[N])
m     net magnetisation
h     hamiltonian
M     magnetisation (float[T])
U     return (float[T])
P     stock price (float[T])

# Lookup
* search for " [##] " for key changeable params
* lattice Nx,Ny     = 24,24
* cutoff K          = 4
* field coupling A  = 20
* prob control B    = 1
* spin coupling J   = 0.12
* daily trades dt   = 100
*/

// [##] settings
#define seed 0
#define T    10000
#define Nx   24
#define Ny   24
#define N    Nx*Ny
#define K    4.

#define A    20
#define B    1
#define J    0.12

int     dt,t,nb,Nb[N][N];
float   m,h,M[T],U[T],P[T],S[N],V[N],D[N][N];
vec     a1,a2,b1,b2,R[N];

/********************************************************/

/* print */

void print_port(){
	/* print portfolios */
	FILE *f = fopen("out/port.csv","w");
	fprintf(f,"index,value\n");
	for(int i=0; i<N; i++)
		fprintf(f,"%d,%f\n",i,V[i]);
	fclose(f);
}

void print_data(){
	/* print time series data */
	int t0=2000; // discard initial non-eqm data
	FILE *f = fopen("out/data.csv","w");
	fprintf(f,"time,price,magnetisation,return\n");
	for(int t=t0; t<T; t++)
		fprintf(f,"%d,%f,%f,%f\n",t-t0,P[t],M[t],U[t]);
	fclose(f);
}

void print_spin(FILE *f){
	/* print spins */
	int n=0;
	for(int i=0; i<Nx; i++){
		for(int j=0; j<Ny; j++){
			if(S[n]==-1) fprintf(f,"0"); // down
			else fprintf(f,"1"); // up
			if(j<Ny-1) fprintf(f,",");
			n++;
		}
		fprintf(f,"\n");
	}
}

void info(int period){
	/* print a line of info to screen */
	if(t%period==0)
		printf(
		"| progress: %5.1f%% | $: %9.6f | mag: %9.6f | ret: %9.6f |\n",
		100.*t/T,P[t-1],M[t-1],U[t-1]);
}

/********************************************************/

/* model */

void init(){
	t=1;

	// lattice basis vectors
	a1=cart2d(Nx,0);
	a2=cart2d(0,Ny);
	b1=cart2d(1./Nx,0);
	b2=cart2d(0,1./Ny);

	int n=0;
	float x;
	/* random initial spin */
	for(int i=0; i<Nx; i++)
		for(int j=0; j<Ny; j++){
			// spin
			x=uniform(0,1);
			if(x<0.5) S[n]=+1; // buy position
			else S[n]=-1; // sell position

			V[n]=0;
			R[n]=cart2d(i,j); // position on lattice
			n++;
		}

	M[0]=mean(S,N);
	P[0]=1; // initial stock price
}

void update(){
	int i;
	float x,p,m,u;

	// [##] a single MC step (transactions happening in a trade day)
	dt=100;
	for(int time=0; time<dt; time++){
		/* MC: update traders' decisions; single trade day */
		/* in each time step, only a single trade */
		m=mean(S,N); // magnetisation

		i=(int)uniform(0,N-1); // sample a trader
		h=0; // hamiltonian
		for(int n=0; n<nb; n++)
			h+=J*S[Nb[i][n]]; // local alignment
		h-=A*S[i]*fabs(m); // reaction to global market atmosphere

		x=uniform(0,1);
		p=1/(1+exp(-2*B*h));
		if(x<p) {S[i]=+1; V[i]+=P[t-1];} // buy
		else {S[i]=-1; V[i]-=P[t-1];} // sell
	}

	/* record current data */
	M[t]=mean(S,N);
	U[t]=M[t]-M[t-1];
	P[t]=P[t-1]*exp(U[t]);

	t++;
}

void iter(){
	while(t<T){
		/* whole time horizon */
		update();
		info(1000); // simulation log to screen
		// [##] (for animation) print whole spin grid
/*
		if(t%500==0){
			char file[sizeof "out/spin_00000.csv"];
			sprintf(file,"out/spin_%05d.csv",t);
			FILE *f=fopen(file,"w");
			print_spin(f);
			fclose(f);
		}
*/
	}
	normalise(U,T);
}

/********************************************************/

/* nbor list */

vec c_to_f(vec c){
	/* Cartesian to fractional */
	vec f;
	f=v(2);
	f.v[0] = mod0(dot(c,b1),1);
	f.v[1] = mod0(dot(c,b2),1);
	return f;
}

vec f_to_c(vec f){
	/* fractional to Cartesian */
	vec c=v(2);
	for(int i=0; i<2; i++)
		c.v[i]=f.v[0]*a1.v[i]+f.v[1]*a2.v[i];
	return c;
}

void make_nbor(){
	/* make nbor list; cutoff=K */
	vec v;
	float d,eps=1e-2;
	for(int i=0; i<N; i++){
		/* get all distances */
		nb=0;
		for(int j=0; j<N; j++){
			// round to avoid floating point err
			v=f_to_c(c_to_f(minus(R[j],R[i]))); // pbc displacement: i->j
			d=roundsf(mag(v),4); // pbc distance
			D[i][j]=d;
			/* pick nbor */
			if(d>0 && d<K+eps){
				Nb[i][nb]=j;
				nb++;
			}
		}
	}
}

/********************************************************/

int main(){
	srand(seed);
	init();
	make_nbor();
	iter();
	print_data();
	return 0;
}
