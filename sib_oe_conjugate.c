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
    
mxArray   *ya;    
mxArray   *Aa;
mxArray   *Ba;
double *y;
double *A;
double *B;
int k;

Aa = mxCreateDoubleMatrix(1+dteta-dnum,1,mxREAL);
A = mxGetPr(Aa);
Ba = mxCreateDoubleMatrix(delay+dnum,1,mxREAL);
B = mxGetPr(Ba);

ya = mxCreateDoubleMatrix(du,1,mxREAL);
y = mxGetPr(ya);

/* Monta Vetores A e B
 */

for (k = 0; k < delay+dnum; k++)
{
    if (k<delay)
    {
        B[k]=0;
    }else{
        B[k]=teta[k-delay];
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

/* Filtra
 */
filt(A,1+dteta-dnum,B,delay+dnum,u,du,y);

J[0]=0;

/* Calcula Custo
 */
for (k = 0; k < du; k++)
{
    J[0]+=(ym[k]-y[k])*(ym[k]-y[k]);
}

J[0]=J[0]/du;

mxDestroyArray(ya);
mxDestroyArray(Aa);
mxDestroyArray(Ba);
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
mxArray   *Aa;
mxArray   *Ba;
mxArray   *ya;    
mxArray   *ea;    
mxArray   *numa;    
mxArray   *dena;    
double *e;
double *A;
double *B;
double *y;
double *num;
double *den;
int j,k;

Aa = mxCreateDoubleMatrix(dteta-dnum+1,1,mxREAL);
A = mxGetPr(Aa);
Ba = mxCreateDoubleMatrix(delay+dnum,1,mxREAL);
B = mxGetPr(Ba);
numa = mxCreateDoubleMatrix(delay+dnum+1,1,mxREAL);
num = mxGetPr(numa);
dena = mxCreateDoubleMatrix((dteta-dnum+1)*2-1,1,mxREAL);
den = mxGetPr(dena);
ya = mxCreateDoubleMatrix(du,1,mxREAL);
y = mxGetPr(ya);
ea = mxCreateDoubleMatrix(du,1,mxREAL);
e = mxGetPr(ea);

/* Calcula o custo
 */

for (k = 0; k < (delay+dnum); k++)
{
    if (k<delay)
    {
        B[k]=0;
    }else{
        B[k]=teta[k-delay];
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
    
filt(A,dteta-dnum+1,B,delay+dnum,u,du,y);
J[0]=0;
    for (k = 0; k < du; k++)
    {
    e[k]=(y[k]-ym[k]);
    J[0]+=(e[k])*(e[k]);
    }
J[0]=J[0]/du;

/* Calcula gradiente com relacao ao numerador
 */

for (j = 0; j < dnum; j++)
{
    for (k = 0; k < (delay+dnum); k++)
    {
        if (k==delay+j)
        {
            num[k]=1;
        }else{
            num[k]=0;;
        }
    }
    for (k = 0; k < dteta-dnum+1; k++)
    {
        if (k==0)
        {
            den[k]=1;
        }else{
            den[k]=teta[k+dnum-1];
        }
    }
    
    filt(den,(dteta-dnum+1),num,(delay+dnum),u,du,y);
    gra[j]=0;
    
    for (k = 0; k < du; k++)
    {
        gra[j]+=(e[k])*(y[k]);
    }
}

for (k = 0; k < (delay+dnum+1); k++)
{
    if (k<delay+1)
    {
        num[k]=0;
    }else{
        num[k]=-teta[k-delay-1];
    }
}

conv(A,dteta-dnum+1,den);
filt(den,(dteta-dnum+1)*2-1,num,delay+dnum+1,u,du,y);

for (j = 0; j < (dteta-dnum); j++)
{
gra[dnum+j]=0;
    for (k = 0; k < (du-j); k++)
    {
        gra[dnum+j]+=(e[k+j])*(y[k]);
    }

}

for (k = 0; k < dteta; k++)
{
    gra[k]=gra[k]/du;
}


mxDestroyArray(Aa);
mxDestroyArray(Ba);
mxDestroyArray(ya);
mxDestroyArray(ea);
mxDestroyArray(numa);
mxDestroyArray(dena);

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
    
    if((i==1)&(j==1))
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


