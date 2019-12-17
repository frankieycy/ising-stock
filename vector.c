#ifndef VECTOR
#define VECTOR

/**** def ****/

typedef struct vec{
	int l;
	float *v;
} vec;

vec v(int l){
	/* init empty vec */
	vec x={l,malloc(sizeof(float)*l)};
	return x;
}

/********************************************************/

/**** handling ****/

vec cart2d(float a, float b){
	/* Cartesian 2d vec */
	vec x=v(2);
	x.v[0]=a; x.v[1]=b;
	return x;
}

/********************************************************/

/**** dynamic arrays ****/

int **int2d(int r, int c){ 
	/* init 2d int array */
	int **x;
	x=malloc(sizeof(int*)*r);
	for(int i=0; i<r; i++) x[i]=malloc(sizeof(int)*c);
	return x;
}

float *float1d(int r){
	/* init 1d float array */
	float *x;
	x=malloc(sizeof(float)*r);
	return x;
}

float **float2d(int r, int c){ 
	/* init 2d float array */
	float **x;
	x=malloc(sizeof(float*)*r);
	for(int i=0; i<r; i++) x[i]=malloc(sizeof(float)*c);
	return x;
}

vec *vec1d(int l){
	/* init 1d vec array */
	vec *x;
	x=malloc(sizeof(vec)*l);
	return x;
}

vec **vec2d(int r, int c){
	/* init 2d vec array */
	vec **x;
	x=malloc(sizeof(vec*)*r);
	for(int i=0; i<r; i++) x[i]=malloc(sizeof(vec)*c);
	return x;
}

/********************************************************/

/**** arithmetic ****/

float mag(vec a){
	/* vector norm */
	float l2=0;
	for(int i=0; i<a.l; i++) l2 += a.v[i]*a.v[i];
	return sqrt(l2);
}

vec add(vec a, vec b){
	/* vector addition */
	if(a.l!=b.l){
		fprintf(stderr,
		"Error in add(vec a, vec b): a and b not of same dimension.\n");
	}
	vec x=v(a.l);
	for(int i=0; i<a.l; i++) x.v[i]=a.v[i]+b.v[i];
	return x;
}

vec minus(vec a, vec b){
	/* vector substraction */
	if(a.l!=b.l){
		fprintf(stderr,
		"Error in minus(vec a, vec b): a and b not of same dimension.\n");
	}
	vec x=v(a.l);
	for(int i=0; i<a.l; i++) x.v[i]=a.v[i]-b.v[i];
	return x;
}

vec vmul(vec a, float b){
	/* scalar multiplication */
	vec x=v(a.l);
	for(int i=0; i<a.l; i++) x.v[i] = a.v[i]*b;
	return x;
}

vec vdiv(vec a, float b){
	/* scalar division */
	return vmul(a,1/b);
}

float dot(vec a, vec b){
	/* dot product */
	if(a.l!=b.l){
		fprintf(stderr,
		"Error in dot(vec a, vec b): a and b dimensions not matched.\n");
		exit(-1);
	}
	float p=0;
	for(int i=0; i<a.l; i++) p += a.v[i]*b.v[i];
	return p;
}

void printv(vec a, char *name){
	printf("%s[] = {",name);
	for(int i=0; i<a.l; i++){
		printf(" %9.4lf ",a.v[i]);
		if(i<a.l-1) printf(",");
	}
	printf("}\n");
}

#endif
