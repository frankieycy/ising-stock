#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "tool.c"

/*
# Model
random walk model of stock

# Definitions
T     total time
t     current time (day)
mu    rate of return (expected daily return)
sig   volatility (std dev of daily return)
U     return (float[T])
P     stock price (float[T])
*/

#define T   10000
#define mu  0
#define sig 0.002

int   t;
float U[T],P[T];

/***********************************************************/

void init(){
	t=1;
	P[0]=1;
}

void update(){
	float z=normal(0,1);
	U[t]=mu+sig*z;
	P[t]=P[t-1]*exp(U[t]);
	t++;
}

void iter(){
	while(t<T) update();
}

void print_data(){
	FILE *f=fopen("out/rw.csv","w");
	int t0=2000;
	for(t=t0; t<T; t++)
		fprintf(f,"%d,%f,%f\n",t-t0,P[t],U[t]);
	fclose(f);
}

/***********************************************************/

int main(){
	init();
	iter();
	print_data();
	return 0;
}
