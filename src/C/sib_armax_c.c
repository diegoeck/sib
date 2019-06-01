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
    double *y0;
    double *y1;
    double *y2;
    double *y3;
    int i,j,k;

    double *pB;
    double *pA;
    double *pC;

    double p1[1]={1};
    double p1n[1]={-1};
    
    pB=malloc( (dB)*sizeof(double) );
    pA=malloc( (dA+1)*sizeof(double) );
    pC=malloc( (dC+1)*sizeof(double) );

    
    y0=malloc((du)*sizeof(double));
    y1=malloc((du)*sizeof(double));
    y2=malloc((du)*sizeof(double));
    y3=malloc((du)*sizeof(double));
    e=malloc((du)*sizeof(double));

    
    for (k = 0; k < dB; k++)
    {
        pB[k]=teta[k];
    }
    
    pA[0] = 1;
    for (k = 0; k < dA; k++)
    {
        pA[k+1]=teta[dB+k];
    }

    pC[0] = 1;
    for (k = 0; k < dC; k++)
    {
        pC[k+1]=teta[dB+dA+k];
    }
    
    filter(pB, pC, dB, dC+1, delay, u, du, y0); 
    filter(pA, pC, dA+1, dC+1, 0, ym, du, y1); 
    
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
        
        filter(p1, pC, 1, dC+1, delay, u, du, y1); 
        for (j = 0; j < dB; j++)
        {
            gra[j]=0;
            for (k = j; k < du; k++)
            {
                gra[j]+=(e[k])*(y1[k-j]);
            }
        }

        // Calcula gradiente com relacao a A
        filter(p1n, pC, 1, dC+1, 1, ym, du, y2); 
        for (j = 0; j < dA; j++)
        {
            gra[dB+j]=0;
            for (k = j; k < du; k++)
            {
                gra[dB+j]+=(e[k])*(y2[k-j]);
            }
        }
        
        // Calcula gradiente com relacao a C
        filter(p1n, pC, 1, dC+1, 1, e, du, y3); 
        for (j = 0; j < dC; j++)
        {
            gra[dA+dB+j]=0;
            for (k = j; k < du; k++)
            {
                gra[dA+dB+j]+=(e[k])*(y3[k-j]);
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
        
         
        
        for (j = 0; j < dB; j++)
        {
            for (k = j; k < dB; k++)
            {
                for (i = 0; i < du; i++)
                {   
                    if ((i>=j)&(i>=k))
                    {
                        H[j][k]+=y1[i-j]*y1[i-k];
                    }
                }
                H[k][j]=H[j][k];
            }
            for (k = 0; k < dA; k++)
            {
                for (i = 0; i < du; i++)
                {   
                    if ((i>=j)&(i>=k))
                    {
                        H[j][k+dB]+=y1[i-j]*y2[i-k];
                    }
                }
                H[k+dB][j]=H[j][k+dB];
            }
            for (k = 0; k < dC; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[j][k+dB+dA]+=y1[i-j]*y3[i-k];
                    }
                }
                H[k+dB+dA][j]=H[j][k+dB+dA];
            }
        }

        for (j = 0; j < dA; j++)
        {
            for (k = j; k < dA; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[j+dB][k+dB]+=y2[i-j]*y2[i-k];
                    }
                }
                H[k+dB][j+dB]=H[j+dB][k+dB];
            }
            for (k = 0; k < dC; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[j+dB][k+dB+dA]+=y2[i-j]*y3[i-k];
                    }
                }
                H[k+dB+dA][j+dB]=H[j+dB][k+dB+dA];
            }
        }
        
        for (j = 0; j < dC; j++)
        {
            for (k = j; k < dC; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[j+dB+dA][k+dB+dA]+=y3[i-j]*y3[i-k];
                    }
                }
                H[k+dB+dA][j+dB+dA]=H[j+dB+dA][k+dB+dA];
            }
        }
           
        

    }
    
    free(y0);
    free(y1);
    free(y2);
    free(y3);
    free(e);

    free(pA);
    free(pB);
    free(pC);

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

    //mexPrintf(" %d %d %d %d  \n",dteta,dA,dB,dC);

    memcpy(tetafim,teta,sizeof(double) *dteta);
    sib_steepest(tetafim);
    //memcpy(teta,tetafim,sizeof(double) *dteta);
    sib_newton(tetafim);
    
    free(gra);
    free(H[0]);
    free(H);
} 

        

        
        
        
        
