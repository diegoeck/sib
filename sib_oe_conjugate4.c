#include <stdint.h>
typedef uint16_t char16_t;
#include "mex.h"
#include <math.h>
#include <string.h>
#include <float.h>
#include <sib_basic.h>
#include <sib_optimize.h>



double *teta;
double *u;
double *ym;
double *tetafim;
int du,dteta,num,delay;




void cost(double* teta, double* J)
{
/* Calcula o Valor Medio Quadratico de (ym(t)-y(t))
 * A(z)y(t)=B(z)u(t)
 * A=[1 teta[dnum] ... teta[dteta-1] ]
 * B=[0 ...(delay)... teta[0] ... teta[dnum-1] 
 */
    
double *y;
int k;
y=malloc((du)*sizeof(double));

// Filtra
filtra(teta,num,dteta-num, u,du,y);

J[0]=0;

// Calcula Custo
for (k = 0; k < du-delay; k++)
{
    J[0]+=(y[k+delay]-y[k])*(y[k+delay]-y[k]);
}

J[0]=J[0]/du;
free(y);
}



void grad(double* teta, int dteta, int delay, int dnum, double* u, int du,double* ym, double* J,double* gra)
{
// Calcula o gradiente da estrutura OE e o custo atual
double *e;
double *A;
double *y;
double *y0;
int j,k;

A=malloc((dteta-dnum+1)*sizeof(double));
y0=malloc((du)*sizeof(double));
y=malloc((du)*sizeof(double));
e=malloc((du)*sizeof(double));


// Calcula o custo
filtra(teta,dnum,dteta-dnum, u,du,y0);

J[0]=0;
for (k = 0; k < du-delay; k++)
{
    e[k]=(y0[k]-ym[k+delay]);
    J[0]+=(e[k])*(e[k]);
}
J[0]=J[0]/du;

// Calcula gradiente com relacao ao numerador

A[0]=1;
for (k = 0; k < dteta-dnum; k++)
{
    A[k+1]=teta[k+dnum];
}
filtra(A,1,dteta-dnum, u,du,y);

for (j = 0; j < dnum; j++)
{
    gra[j]=0;
    for (k = 0; k < du-delay-j; k++)
    {
        gra[j]+=(e[k+j])*(y[k]);
    }
}

// Calcula gradiente com relacao ao denominador
filtra(A,1,dteta-dnum, y0,du,y);
for (j = 0; j < (dteta-dnum); j++)
{
    gra[j+dnum]=0;
    for (k = 0; k < du-delay-j-1; k++)
    {
        gra[j+dnum]-=(e[k+j+1])*(y[k]);
    }
}


for (k = 0; k < dteta; k++)
{
    gra[k]=gra[k]/du;
}


free(y);
free(y0);
free(e);
free(A);

}




void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{

du = mxGetM(prhs[0]);
dteta = mxGetM(prhs[2]);

u = mxGetPr(prhs[0]);
ym = mxGetPr(prhs[1]);
teta = mxGetPr(prhs[2]);
delay = *mxGetPr(prhs[3]);
num = *mxGetPr(prhs[4]);

plhs[0] = mxCreateDoubleMatrix(dteta,1,mxREAL);
tetafim = mxGetPr(plhs[0]);

sib_steepest(tetafim);

} 


