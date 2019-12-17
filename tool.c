#ifndef TOOL
#define TOOL

/**** misc ****/

void line(){
	printf("-------------------------------------\n");
}

void fline(FILE *f){
	fprintf(f,"-------------------------------------\n");
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

/**** handling ****/

float window(float *a, int n, int m, float func(float*,int)){
	/* operate func on a[n] to a[m] */
	return func(a+n,m-n+1);
}

/********************************************************/

/**** statistics ****/

float sum(float *a, int l){
	float x=0;
	for(int i=0; i<l; i++) x+=a[i];
	return x;
}

float mean(float *a, int l){
	return sum(a,l)/l;
}

float wmean(float *a, int l){
	/* linearly weighted mean */
	float x=0;
	for(int i=0; i<l; i++) x+=(i+1)*a[i];
	return x/(l*(l+1)/2);
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
	/* convert to mean 0, variance 1 */
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

/**** tech indicators ****/

float SMA(float *a, int l){
	/* simple moving average */
	return mean(a,l);
}

float WMA(float *a, int l){
	/* weighted moving average */
	return wmean(a,l);
}

float RSI(float *a, int l){
	/* relative strength index */
	float x,R=0,D=0;
	for(int i=1; i<l; i++){
		x=a[i]-a[i-1];
		if(x>0) R+=x; // rise
		if(x<0) D-=x; // drop
	}
	if(D==0) return 100;
	return 100*(1-1/(1+R/D));
}

/********************************************************/

/**** sampling ****/

float uniform(float a, float b){
	float u=(float)rand()/RAND_MAX;
	return a+(b-a)*u;
}

float normal(float a, float b){
	float u1=uniform(0,1),u2=uniform(0,1),
	z=sqrt(-2*log(u1))*cos(2*M_PI*u2);
	return a+b*z;
}

int poisson(float a){
	int k=0;
	float u,L=exp(-a),p=1;
	while(p>L){
		u=uniform(0,1);
		p*=u;
		k++;
	}
	return k-1;
}

#endif
