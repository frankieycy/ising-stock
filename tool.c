#ifndef TOOL
#define TOOL
#include <stdlib.h>
#include <math.h>

/* misc */

void line(){
	printf("--------------\n");
}

void fline(FILE *f){
	fprintf(f,"--------------\n");
}

float mod0(float x, float d){
	/* modulus mapped to [-d/2,d/2] */
	while(x>+d/2) x-=d;
	while(x<-d/2) x+=d;
	return x;
}

float min0(float a, float b){
	/* min of two num */
	if(a<b) return a;
	return b;
}

float roundsf(float a, int n){
	/* round to n sig fig */
	return roundf(a*pow(10,n))/pow(10,n);
}

/********************************************************/

/* statistics */

float sum(float *a, int l){
	float x=0;
	for(int i=0; i<l; i++) x+=a[i];
	return x;
}

float mean(float *a, int l){
	return sum(a,l)/l;
}

float var(float *a, int l){
	/* variance */
	float x=0;
	for(int i=0; i<l; i++) x+=a[i]*a[i];
	return x/l-pow(mean(a,l),2);
}

float std(float *a, int l){
	/* standard deviation */
	return sqrt(var(a,l));
}

void normalise(float *a, int l){
	float m,s;
	m=mean(a,l);
	s=std(a,l);
	for(int i=0; i<l; i++) a[i]=(a[i]-m)/s;
}

float min(float *a, int l){
	float x=a[0];
	for(int i=1; i<l; i++) if(a[i]<x) x=a[i];
	return x;
}

float max(float *a, int l){
	float x=a[0];
	for(int i=1; i<l; i++) if(a[i]>x) x=a[i];
	return x;
}

void print_stat(float *a, int l){
	/* summary statistics */
	printf("min  : %.4f\n"
	       "max  : %.4f\n"
	       "mean : %.4f\n"
	       "std  : %.4f\n"
	,min(a,l),max(a,l),mean(a,l),std(a,l));
}

/********************************************************/

/* sampling */

float uniform(float a, float b){
	float u=(float)rand()/RAND_MAX;
	return a+(b-a)*u;
}

float normal(float a, float b){
	float u1=uniform(0,1),u2=uniform(0,1),
	z=sqrt(-2*log(u1))*cos(2*M_PI*u2);
	return a+b*z;
}

#endif
