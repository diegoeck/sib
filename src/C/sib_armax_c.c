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
    int i,j,k;

    double *e;
    double *yf;
    double *dydA;
    double *dydB;
    double *dydC;

    double *pA;
    double *pB;
    double *pC;

    double p1[1]={1};
    double p1n[1]={-1};
    
    pA = malloc( (dA+1)*sizeof(double) );
    pB = malloc( (dB)*sizeof(double) );
    pC = malloc( (dC+1)*sizeof(double) );
    
    yf = malloc((du)*sizeof(double));
    dydA = malloc((du)*sizeof(double));
    dydB = malloc((du)*sizeof(double));
    dydC = malloc((du)*sizeof(double));
    e = malloc((du)*sizeof(double));

    pA[0] = 1;
    for (k = 0; k < dA; k++)
    {
        pA[k+1] = teta[k];
    }

    for (k = 0; k < dB; k++)
    {
        pB[k] = teta[dA+k];
    }
    

    pC[0] = 1;
    for (k = 0; k < dC; k++)
    {
        pC[k+1] = teta[dA+dB+k];
    }
    
    filter(pB, pC, dB, dC+1, delay, u, du, e); 
    filter(pA, pC, dA+1, dC+1, 0, ym, du, yf); 
    
    J[0]=0;
    for (k = 0; k < du; k++)
    {
        e[k] = (e[k]-yf[k]);
        J[0] += (e[k])*(e[k]);
    }
    J[0] = sqrt(J[0]/du);

    
    
    if(mode==1 | mode==2){
        
        
        for (j = 0; j < dA+dB+dC; j++)
        {
            gra[j] = 0;
        }

        // Calcula gradiente com relacao a A
        filter(p1n, pC, 1, dC+1, 1, ym, du, dydA); 
        for (j = 0; j < dA; j++)
        {
            for (k = j; k < du; k++)
            {
                gra[j] += (e[k])*(dydA[k-j]);
            }
        }
        
        
        // Calcula gradiente com relacao a B
        filter(p1, pC, 1, dC+1, delay, u, du, dydB); 
        for (j = 0; j < dB; j++)
        {
            for (k = j; k < du; k++)
            {
                gra[dA+j] += (e[k])*(dydB[k-j]);
            }
        }

        
        // Calcula gradiente com relacao a C
        filter(p1n, pC, 1, dC+1, 1, e, du, dydC); 
        for (j = 0; j < dC; j++)
        {
            for (k = j; k < du; k++)
            {
                gra[dA+dB+j] += (e[k])*(dydC[k-j]);
            }
        }
        
    }
    
    if(mode==2)
    {
    
        for (j = 0; j < dteta; j++)
        {
            for (k = 0; k < dteta; k++)
            {
                H[j][k] = 0;
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
                        H[j][k] += dydA[i-j]*dydA[i-k];
                    }
                }
                H[k][j] = H[j][k];
            }
            for (k = 0; k < dB; k++)
            {
                for (i = 0; i < du; i++)
                {   
                    if ((i>=j)&(i>=k))
                    {
                        H[j][dA+k] += dydA[i-j]*dydB[i-k];
                    }
                }
                H[k+dA][j] = H[j][k+dA];
            }
            for (k = 0; k < dC; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[j][dA+dB+k] += dydA[i-j]*dydC[i-k];
                    }
                }
                H[dA+dB+k][j] = H[j][dA+dB+k];
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
                        H[dA+j][dA+k] += dydB[i-j]*dydB[i-k];
                    }
                }
                H[dA+k][dA+j] = H[dA+j][dA+k];
            }
            for (k = 0; k < dC; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[dA+j][dA+dB+k] += dydB[i-j]*dydC[i-k];
                    }
                }
                H[dA+dB+k][dA+j] = H[dA+j][dA+dB+k];
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
                        H[dA+dB+j][dA+dB+k] += dydC[i-j]*dydC[i-k];
                    }
                }
                H[dA+dB+k][dA+dB+j] = H[dA+dB+j][dA+dB+k];
            }
        }
           
        

    }
    
    free(e);
    free(yf);
    free(dydA);
    free(dydB);
    free(dydC);

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
    
    if (dteta == dA+dB+dC)
    {
        gra=malloc((dteta)*sizeof(double));
        H = malloc(dteta * sizeof(double *));
        H[0] = malloc(dteta * dteta * sizeof(double));
        for(i = 1; i < dteta; i++)
            H[i] = H[0] + i * dteta;

        plhs[0] = mxCreateDoubleMatrix(dteta,1,mxREAL);
        tetafim = mxGetPr(plhs[0]);

        memcpy(tetafim,teta,sizeof(double) *dteta);
        sib_steepest(tetafim);
        sib_newton(tetafim);
    
        free(gra);
        free(H[0]);
        free(H);
    }else{
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        tetafim = mxGetPr(plhs[0]);
        tetafim = 0;
        
        mexPrintf("ERROR: Length of theta should be na+nb+nc !!! \n");
    }
} 

        

        
        
        
        
