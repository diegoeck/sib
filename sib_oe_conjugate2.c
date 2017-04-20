#include <stdint.h>
typedef uint16_t char16_t;
#include "mex.h"
#include <math.h>
#include <string.h>
#include <float.h>

void filt(double* A,int dA, double* B,int dB, double* u, int du, double* y)
{
/* Funcao Filtro A(z)y(t)=B(z)u(t)
 *
 * Considera que A[0]=1 
 */
    
int i;
int j;
int k;
for (i = 0; i <du; i++)  
{
    y[i]=0;
    for (j = 0; j < dB; j++)
    {
        if ((i-j)>=0)
        {
          y[i] += u[i-j]*B[j]; 
        }
    }
    for (k = 1; k < dA; k++)
    {
        if ((i-k)>=0)
        {
          y[i] -= y[i-k]*A[k]; 
        }
    }
}

}

void filtra(double* teta,int n,int d, double* u, int du, double* y)
{
/* Funcao Filtro A(z)y(t)=B(z)u(t)
 *
 * Considera que A[0]=1 
 */
    
int i;
int j;
for (i = 0; i <du; i++)  
{
    y[i]=0;
    for (j = 0; j < n; j++)
    {
        if ((i-j)>=0)
        {
          y[i] += u[i-j]*teta[j]; 
        }
    }
    for (j = 0; j < d; j++)
    {
        if ((i-j-1)>=0)
        {
          y[i] -= y[i-j-1]*teta[j+n]; 
        }
    }
}

}



double mult_scalar(double* A, double* B, int tamanho)
{
/* Multiplicacao escalar de vetores
 */

double resultado;    
int i;
resultado=0;
for (i = 0; i <tamanho; i++)  
{
    resultado+=A[i] * B[i];
}
return resultado;
}

void cost(double* teta, int dteta, int delay, int dnum, double* u, int du,double* ym, double* J)
{
/* Calcula o Valor Medio Quadratico de (ym(t)-y(t))
 * A(z)y(t)=B(z)u(t)
 * A=[1 teta[dnum] ... teta[dteta-1] ]
 * B=[0 ...(delay)... teta[0] ... teta[dnum-1] 
 */
    
double *y;
int k;
y=malloc((du)*sizeof(double));


/* Monta Vetores A e B
 */

/* Filtra
 */
filtra(teta,dnum,dteta-dnum, u,du,y);

J[0]=0;

/* Calcula Custo
 */
for (k = 0; k < du; k++)
{
    J[0]+=(ym[k]-y[k])*(ym[k]-y[k]);
}

J[0]=J[0]/du;
free(y);
}

void conv(double* A,int dA, double* B)
{
/* Convolucao de A com A
 * Retorna em B
 */
int i;
int j;

for (i = 0; i < (dA*2-1); i++)
{
    B[i]=0;
}

for (i = 0; i < dA; i++)
{
    for (j = 0; j < dA; j++)
    {
        B[i+j]+=A[i]*A[j];
    }
}

}

void subtr(double* A,double* B, int dA, double* C)
{
/* Subtracao vetorial
  */
int i;

for (i = 0; i < dA; i++)
{
    C[i]=A[i]-B[i];
}

}

double norm(double* A,int dA)
{
/* Norma de A
 */
return sqrt(mult_scalar(A,A,dA));
}

void grad(double* teta, int dteta, int delay, int dnum, double* u, int du,double* ym, double* J,double* gra)
{
/* Calcula o gradiente da estrutura OE e o custo atual
 */
double *e;
double *A;
double *B;
double *y;
double *y0;
double *num;
double *den;
int j,k;


num=malloc((delay+dnum+1)*sizeof(double));
den=malloc(((dteta-dnum+1)*2-1)*sizeof(double));
A=malloc((dteta-dnum+1)*sizeof(double));
y0=malloc((du)*sizeof(double));
y=malloc((du)*sizeof(double));
e=malloc((du)*sizeof(double));


/* Calcula o custo
 */

filtra(teta,dnum,dteta-dnum, u,du,y0);

    
//filt(A,dteta-dnum+1,B,delay+dnum,u,du,y0);
J[0]=0;
for (k = 0; k < du; k++)
{
    e[k]=(y0[k]-ym[k]);
    J[0]+=(e[k])*(e[k]);
}
J[0]=J[0]/du;

/* Calcula gradiente com relacao ao numerador
 */

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
        gra[j]+=(e[k+delay+j])*(y[k]);
    }
}
/*
for (k = 0; k < (delay+dnum+1); k++)
{
    if (k<delay+1)
    {
        num[k]=0;
    }else{
        num[k]=-teta[k-delay-1];
    }
}
for (k = 0; k < dteta-dnum+1; k++)
{
    if (k==0)
    {
        A[k]=1;
    }else{
        A[k]=teta[k+dnum-1];
    }
}

conv(A,dteta-dnum+1,den);
filt(den,(dteta-dnum+1)*2-1,num,delay+dnum+1,u,du,y);
*/
filtra(A,1,dteta-dnum, y0,du,y);
for (j = 0; j < (dteta-dnum); j++)
{
    gra[j+dnum]=0;
    for (k = 0; k < du-delay-j-1; k++)
    {
        gra[j+dnum]-=(e[k+delay+j+1])*(y[k]);
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
free(num);
free(den);

}

void oe_steepest(double* tetainicio, int dteta, int delay, int tnum, double* u, int du,double* ym,double* tetafim)
{
mxArray   *tetatempa;    
mxArray   *graa;
mxArray   *dira;

int i,j,k;
double *tetatemp;
double *gra;
double *dir;
double passo;
double N;
double J1[1];
double J2[1];



tetatempa = mxCreateDoubleMatrix(dteta,1,mxREAL);
tetatemp = mxGetPr(tetatempa);
graa = mxCreateDoubleMatrix(dteta,1,mxREAL);
gra = mxGetPr(graa);
dira = mxCreateDoubleMatrix(dteta,1,mxREAL);
dir = mxGetPr(dira);

memcpy(tetafim,tetainicio,sizeof(double) *dteta);

    
    passo=0.00001;


for(i=0;i<1000;i++)
{

for(j=0;j<1000;j++)
{
    
    grad(tetafim,dteta,delay,tnum,u,du,ym,J1,gra);
    
    if((i==0)&(j==0))
    {
    memcpy(dir,gra,sizeof(double) *dteta);
    }

    for (k = 0; k < dteta; k++)
    {
        dir[k]=(4*dir[k]+gra[k])/5;
    }
    
    N=norm(dir,dteta); 

    for (k = 0; k < dteta; k++)
    {
        tetatemp[k]=tetafim[k]-passo*dir[k]/N;
    }
    
    cost(tetatemp,dteta,delay,tnum, u, du, ym, J2);
        
        if (J2[0]>J1[0])
        {
            passo=passo*0.99;

        }else{   
            passo=passo*1.01;
            memcpy(tetafim,tetatemp,sizeof(double) *dteta);

        }
    
    
}    
    mexPrintf("%d %1.10f %1.10f\n",i,J1[0],passo);

    if (passo<0.0000001)
    {
        break;
    }    
    
    

}



mxDestroyArray(tetatempa);
mxDestroyArray(graa);

}
 




void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
double *teta;
double *u;
double *ym;
double *tetafim;
double *delay;
double *num;
double N[1];
double tt[16];
int du,dteta;

du = mxGetM(prhs[0]);
dteta = mxGetM(prhs[2]);

u = mxGetPr(prhs[0]);
ym = mxGetPr(prhs[1]);
teta = mxGetPr(prhs[2]);
delay = mxGetPr(prhs[3]);
num = mxGetPr(prhs[4]);


plhs[0] = mxCreateDoubleMatrix(dteta,1,mxREAL);
tetafim = mxGetPr(plhs[0]);


oe_steepest(teta,dteta,(int)*delay,(int)*num,u,du,ym,tetafim);

} 


