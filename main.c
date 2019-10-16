#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include "tool.c"
#include "vector.c"

/*
Note:
	- Last: 17/10/2019

Folders:
	- `out`: data outputs
	- `plt`: python plots
	- `bak`: backups
	- `rw`: random walk model (later)

Model:
	- Ising stock market
	- market: fundamentalists, interacting traders
	- interacting traders placed in 2d lattice (2d Ising)
	- spin: buy/sell decision
	- locally seeking ferromagnetic order (nbor-interaction)
	- indiv interacting traders influenced by nbor
	- globally escaping ferromagnetic order (external field)
	- returns & stock prices updated by net magnetisation

Definitions:
	T:		total time
	dt:		num of MC sim between each time step
	t:		current time

	N:		num of interacting traders
	K:		cutoff distance (to pick nbor)
	R:		lattice coord (vec[N][N])
	nb:		num of nbor
	Nb:		nbor list (int[N][N])
	D:		pbc nbor distances (float[N][N])

	A:		external field strength
	B:		controls state transition
	J:		coupling strength

	S:		spin (+1 or -1) (float[N])
	m:		net magnetisation
	h:		hamiltonian
	H:		hamiltonian (float[T])
	M:		magnetisation (float[T])
	U:		return (float[T])
	P:		stock price (float[T])

	f0:		data
	f1:		spin
*/

#define T  10000
#define dt 100
#define Nx 32
#define Ny 32
#define K  1.

#define A 20
#define B 1
#define J 1

int   t,N,nb,**Nb;
float m,h,H[T],M[T],U[T],P[T],*S,**D;
vec   a1,a2,b1,b2,*R;
FILE  *f0,*f1;

/********************************************************/

void update(){
	int i;
	float x,p,m,u;

	for(int n=0; n<dt; n++){
		/* MC: update traders' decisions */
		m=mean(S,N);

		i=(int)uniform(0,N); // sample a trader
		h=0; // hamiltonian
		for(int j=0; j<nb; j++)
			h+=J*S[Nb[i][j]];
		h-=A*S[i]*fabs(m);

		x=uniform(0,1);
		p=1/(1+exp(-2*B*h));
		if(x<p) S[i]=+1; // buy
		else S[i]=-1; // sell
	}

	/* record data */
	H[t]=h;
	M[t]=m;
	U[t]=M[t]-M[t-1];
	P[t]=P[t-1]*exp(U[t]);

	t++;
}

void print_data(FILE *f){
	/* print to data file */
	int t0=2000; // discard initial data
	fprintf(f,"time,price,magnetisation,return\n");
	for(int t=t0; t<T; t++)
		fprintf(f,"%d,%f,%f,%f\n",t-t0,P[t],M[t],U[t]);
}

void print_spin(FILE *f){
	/* print to spin file */
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

void iter(){
	while(t<T){
		update();
		if(t%1000==0){
			char file[sizeof "out/spin_00000.csv"];
			sprintf(file,"out/spin_%05d.csv",t);
			f1=fopen(file,"w");
			print_spin(f1);
			fclose(f1);
		}
	}
	normalise(U,T);
}

/********************************************************/

void init(){
	t=1;
	N=Nx*Ny;
	M[0]=0;
	P[0]=1; // initial stock price

	Nb=int2d(N,N);
	S=float1d(N);
	D=float2d(N,N);
	R=vec1d(N);

	// lattice vectors
	a1=cart2d(Nx,0);
	a2=cart2d(0,Ny);
	b1=cart2d(1./Nx,0);
	b2=cart2d(0,1./Ny);

	int n=0;
	float x;
	/* random initial spin */
	for(int i=0; i<Nx; i++)
		for(int j=0; j<Ny; j++){
			x=uniform(0,1);
			if(x<0.5) S[n]=+1;
			else S[n]=-1;
			R[n]=cart2d(i,j); // position on lattice
			n++;
		}
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
	srand(0);
	init();
	make_nbor();
	iter();

	f0=fopen("out/data.csv","w");
	print_data(f0);
	fclose(f0);
	return 0;
}
