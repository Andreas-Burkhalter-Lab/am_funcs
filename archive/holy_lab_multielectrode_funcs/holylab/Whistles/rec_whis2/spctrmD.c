// Modified version of NR's spctrm, Tim Holy 10/13/99.
// Works in memory rather than with a file,
// and output is a double rather than a float
// Useful for interfacing with MATLAB in an API.
// If overlap = 1, input data must contain (2*k+1)*m points
// If overlap = 0, input data must contain 4*k*m points
#include <math.h>
#include <stdio.h>
#define NRANSI
#include "nrutil.h"
#define WINDOW(j,a,b) (1.0-fabs((((j)-1)-(a))*(b)))       /* Bartlett */

extern "C" void spctrmD(float d[], double p[], int m, int k, int ovrlap);
extern "C" void four1(float d[], unsigned long nn, int isign);

void spctrmD(float d[], double p[], int m, int k, int ovrlap)
{
	void four1(float d[], unsigned long nn, int isign);
	int mm,m44,m43,m4,kk,joffn,joff,j2,j;
	long i;
	float w,facp,facm,*w1,*w2,sumw=0.0,den=0.0;

	mm=m+m;
	m43=(m4=mm+mm)+3;
	m44=m43+1;
	w1=vector(1,m4);
	w2=vector(1,m);
	facm=m;
	facp=1.0/m;
	for (j=1;j<=mm;j++) sumw += SQR(WINDOW(j,facm,facp));
	for (j=1;j<=m;j++) p[j]=0.0;
	i = 1;
	if (ovrlap)
		//for (j=1;j<=m;j++) fscanf(fp,"%f",&w2[j]);
		for (j=1;j<=m;j++) w2[j] = d[i++];
	for (kk=1;kk<=k;kk++) {
		for (joff = -1;joff<=0;joff++) {
			if (ovrlap) {
				for (j=1;j<=m;j++) w1[joff+j+j]=w2[j];
				//for (j=1;j<=m;j++) fscanf(fp,"%f",&w2[j]);
				for (j = 1; j <= m; j++) w2[j] = d[i++];
				joffn=joff+mm;
				for (j=1;j<=m;j++) w1[joffn+j+j]=w2[j];
			} else {
				for (j=joff+2;j<=m4;j+=2)
					//fscanf(fp,"%f",&w1[j]);
					w1[j] = d[i++];
			}
		}
		for (j=1;j<=mm;j++) {
			j2=j+j;
			w=WINDOW(j,facm,facp);
			w1[j2] *= w;
			w1[j2-1] *= w;
		}
		four1(w1,mm,1);
		p[1] += (SQR(w1[1])+SQR(w1[2]));
		for (j=2;j<=m;j++) {
			j2=j+j;
			p[j] += (SQR(w1[j2])+SQR(w1[j2-1])
				+SQR(w1[m44-j2])+SQR(w1[m43-j2]));
		}
		den += sumw;
	}
	den *= m4;
	for (j=1;j<=m;j++) p[j] /= den;
	free_vector(w2,1,m);
	free_vector(w1,1,m4);
}
#undef WINDOW
#undef NRANSI
