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

int du;
int dteta;
int dA;
int dB;
int dC;
int delay;
int mode;

void grad(double* teta, double* J)
{
    // Calcula o gradiente da estrutura OE e o custo atual
    double *e;
    double *A;
    double *y;
    double *y0;
    double *y1;
    double *y2;
    double *y3;
    int i,j,k;

    y0=malloc((du)*sizeof(double));
    y1=malloc((du)*sizeof(double));
    y2=malloc((du)*sizeof(double));
    y3=malloc((du)*sizeof(double));
    e=malloc((du)*sizeof(double));

    
    //mexPrintf(" %f %f %f \n",teta[0],teta[1],teta[2]);
    
    A=malloc((dB+dC)*sizeof(double));
    for (k = 0; k < dB; k++)
    {
        A[k]=teta[k];
    }
    for (k = 0; k < dC; k++)
    {
        A[dB+k]=teta[dA+dB+k];
    }
    
    

    
    filtra(A, dB, dC, delay, u, du, y0); 
    free(A);
    
    
    A=malloc((dA+dC+1)*sizeof(double));
    A[0]=1;
    for (k = 0; k < dA+dC+1; k++)
    {
        A[k+1]=teta[dB+k];
    }
    filtra(A, dA+1, dC, 0, ym, du, y1);    
    free(A);
    
    J[0]=0;
    for (k = 0; k < du; k++)
    {
        e[k]=(y0[k]-y1[k]);
        J[0]+=(e[k])*(e[k]);
    }

    
    
    if(mode==1 | mode==2){
        
        
        for (j = 0; j < dA+dB+dC; j++)
        {
            gra[j]=0;
        }
        
        // Calcula gradiente com relacao a B
        A=malloc((1+dC)*sizeof(double));
        A[0]=1;
        for (k = 0; k < dC; k++)
        {
            A[k+1]=teta[dA+dB+k];
        }
        filtra(A, 1, dC, delay, u, du, y1);
        

        for (j = 0; j < dB; j++)
        {
            gra[j]=0;
            for (k = j; k < du; k++)
            {
                gra[j]+=(e[k])*(y1[k-j]);
            }
        }

        // Calcula gradiente com relacao a A
        A[0]=-1;
        filtra(A, 1, dC, 1, ym, du, y2);
        for (j = 0; j < dA; j++)
        {
            gra[dB+j]=0;
            for (k = j; k < du; k++)
            {
                gra[dB+j]+=(e[k])*(y2[k-j]);
            }
        }
        
        
        // Calcula gradiente com relacao a C
        filtra(A, 1, dC, 1, e, du, y3);
        for (j = 0; j < dC; j++)
        {
            gra[dA+dB+j]=0;
            for (k = j; k < du; k++)
            {
                gra[dA+dB+j]+=(e[k])*(y3[k-j]);
            }
        }
        
        //mexPrintf(" %f %f %f \n",gra[0],gra[1],gra[2]);

        
        free(A);

    }
    
    if(mode==2)
    {

        for (j = 0; j < dteta; j++)
        {
            for (k = j; k < dteta; k++)
            {
                H[j][k]=0;
            }
        }
        
        for (i = 0; i < du; i++)
        {            
            for (j = 0; j < dB; j++)
            {
                for (k = j; k < dB; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                        H[j][k]+=(y1[i-k])*(y1[i-j]);
                        H[k][j]=H[j][k];
                    }
                }
                for (k = dB; k < dB+dA; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                        H[j][k]+=(y1[i-j])*(y2[i-k]);
                        H[k][j]=H[j][k];
                    }
                }
                for (k = dB+dA; k < dB+dA+dC; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                        H[j][k]+=(y1[i-j])*(y3[i-k]);
                        H[k][j]=H[j][k];
                    }
                }

            }
        
        
            for (j = dB; j < dB+dA; j++)
            {
                for (k = dB; k < dB+dA; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                     H[j][k]+=(y2[i-j])*(y2[i-k]);
                     H[k][j]=H[j][k];
                    }
                }
                for (k = dB+dA; k < dB+dA+dC; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                     H[j][k]+=(y2[i-j])*(y3[i-k]);
                     H[k][j]=H[j][k];
                    }
                }

            }
    
            
            for (j = dB+dA; j < dB+dA+dC; j++)
            {
                for (k = dB+dA; k < dB+dA+dC; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                     H[j][k]+=(y3[i-j])*(y3[i-k]);
                     H[k][j]=H[j][k];
                    }
                }

            }
    
        
        }

    }
    
    free(y);
    free(y2);
    free(y0);
    free(e);

}




void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{

    int i;
    
    du = mxGetM(prhs[0]);
    dteta = mxGetM(prhs[2]);

    u = mxGetPr(prhs[0]);
    ym = mxGetPr(prhs[1]);
    teta = mxGetPr(prhs[2]);
    dA = *mxGetPr(prhs[3]);
    dB = *mxGetPr(prhs[4]);
    dC = *mxGetPr(prhs[5]);
    delay = *mxGetPr(prhs[6]);
    
    

    gra=malloc((dteta)*sizeof(double));

    H = malloc(dteta * sizeof(double *));
    H[0] = malloc(dteta * dteta * sizeof(double));
    for(i = 1; i < dteta; i++)
        H[i] = H[0] + i * dteta;


    plhs[0] = mxCreateDoubleMatrix(dteta,1,mxREAL);
    tetafim = mxGetPr(plhs[0]);

    mexPrintf(" %d %d %d %d  \n",dteta,dA,dB,dC);

    
    sib_steepest(tetafim);
    memcpy(teta,tetafim,sizeof(double) *dteta);
    sib_newton(tetafim);
} 

        

        
        
        
        
