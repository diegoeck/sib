#include <stdint.h>
#include "mex.h"
#include <math.h>
#include <string.h>
#include <float.h>
typedef uint16_t char16_t;
#include <sib_basic.h>
#include <sib_optimize.h>

double *teta;
double *u;
double *ym;
double *tetafim;
double **H;
double *gra;

int du,dteta,dnum,delay;
int mode;

void grad(double* teta, double* J)
{
    // Calcula o gradiente da estrutura OE e o custo atual
    double *e;
    double *A;
    double *y;
    double *y2;
    double *y0;
    int i,j,k;

    A=malloc((dteta-dnum+1)*sizeof(double));
    y0=malloc((du)*sizeof(double));
    y=malloc((du)*sizeof(double));
    y2=malloc((du)*sizeof(double));
    e=malloc((du)*sizeof(double));


    
    // Calcula o custo
    filtra(teta, dnum, dteta-dnum, delay, u, du, y0);

    J[0]=0;
    for (k = 0; k < du; k++)
    {
        e[k]=(y0[k]-ym[k]);
        J[0]+=(e[k])*(e[k]);
    }

    if(mode==1 | mode==2){
        // Calcula gradiente com relacao ao numerador

        A[0]=1;
        for (k = 0; k < dteta-dnum; k++)
        {
            A[k+1]=teta[k+dnum];
        }
        filtra(A, 1, dteta-dnum, delay, u, du, y);

        for (j = 0; j < dnum; j++)
        {
            gra[j]=0;
            for (k = j; k < du; k++)
            {
                gra[j]+=(e[k])*(y[k-j]);
            }
        }

        // Calcula gradiente com relacao ao denominador
        A[0]=-1;
        filtra(A, 1, dteta-dnum, 1, y0, du, y2);
        for (j = 0; j < (dteta-dnum); j++)
        {
            gra[j+dnum]=0;
            for (k = j; k < du; k++)
            {
                gra[j+dnum]+=(e[k])*(y2[k-j]);
            }
        }
    }
    
    if(mode==2)
    {

        for (j = 0; j < dteta; j++)
        {
            for (k = 0; k < dteta; k++)
            {
                H[j][k]=0;
            }
        }
        
         
        
        for (j = 0; j < dnum; j++)
        {
            for (k = j; k < dnum; k++)
            {
                for (i = 0; i < du; i++)
                {   
                    if ((i>=j)&(i>=k))
                    {
                        H[j][k]+=y[i-j]*y[i-k];
                    }
                }
                H[k][j]=H[j][k];
            }
            for (k = 0; k < dteta-dnum; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[j][k+dnum]+=y[i-j]*y2[i-k];
                    }
                }
                H[k+dnum][j]=H[j][k+dnum];
            }
        }

        for (j = 0; j < dteta-dnum; j++)
        {
            for (k = j; k < dteta-dnum; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[j+dnum][k+dnum]+=y2[i-j]*y2[i-k];
                    }
                }
                H[k+dnum][j+dnum]=H[j+dnum][k+dnum];
            }
        }
        

    }
    
    free(y);
    free(y2);
    free(y0);
    free(e);
    free(A);

}




void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{

    int i;
    
    du = mxGetM(prhs[0]);
    dteta = mxGetM(prhs[2]);

    u = mxGetPr(prhs[0]);
    ym = mxGetPr(prhs[1]);
    teta = mxGetPr(prhs[2]);
    dnum = *mxGetPr(prhs[3]);
    delay = *mxGetPr(prhs[4]);
    
    

    gra=malloc((dteta)*sizeof(double));

    H = malloc(dteta * sizeof(double *));
    H[0] = malloc(dteta * dteta * sizeof(double));
    for(i = 1; i < dteta; i++)
        H[i] = H[0] + i * dteta;


    plhs[0] = mxCreateDoubleMatrix(dteta,1,mxREAL);
    tetafim = mxGetPr(plhs[0]);

    memcpy(tetafim,teta,sizeof(double) *dteta);

    sib_steepest(tetafim);
    //memcpy(teta,tetafim,sizeof(double) *dteta);
    sib_newton(tetafim);
    
    free(gra);
    free(H[0]);
    free(H);

} 

        

        
        
        
        
