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
int dB;
int dF;
int delay;
int mode;

void grad(double* teta, double* J)
{
    // Calcula o gradiente da estrutura OE e o custo atual
    double *e;
    double *y0;
    double *dydB;
    double *dydF;
    int i,j,k;

    double *pB;
    double *pF;

    double p1[1]={1};
    double p1n[1]={-1};
    
    pB = malloc( (dB)*sizeof(double) );
    pF = malloc( (dF+1)*sizeof(double) );

    e=malloc((du)*sizeof(double));
    y0=malloc((du)*sizeof(double));
    dydB=malloc((du)*sizeof(double));
    dydF=malloc((du)*sizeof(double));


    for (k = 0; k < dB; k++)
    {
        pB[k] = teta[k];
    }

    pF[0] = 1;
    for (k = 0; k < dF; k++)
    {
        pF[k+1] = teta[dB+k];
    }
    
    // Calcula o custo
    filter(pB, pF, dB, dF+1, delay, u, du, y0); 
    J[0]=0;
    for (k = 0; k < du; k++)
    {
        e[k]=(y0[k]-ym[k]);
        J[0]+=(e[k])*(e[k]);
    }
    J[0] = sqrt(J[0]/du);


    if(mode==1 | mode==2){
        // Calcula gradiente com relacao ao numerador

        filter(p1, pF, 1, dF+1, delay, u, du, dydB); 
        for (j = 0; j < dB; j++)
        {
            gra[j]=0;
            for (k = j; k < du; k++)
            {
                gra[j]+=(e[k])*(dydB[k-j]);
            }
        }

        // Calcula gradiente com relacao ao denominador
        filter(p1n, pF, 1, dF+1, 1, y0, du, dydF); 
        for (j = 0; j < dF; j++)
        {
            gra[dB+j]=0;
            for (k = j; k < du; k++)
            {
                gra[dB+j] += (e[k])*(dydF[k-j]);
            }
        }
    }
    
    if(mode==2)
    {

        for (j = 0; j < dB+dF; j++)
        {
            for (k = 0; k < dB+dF; k++)
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
                        H[j][k] += dydB[i-j]*dydB[i-k];
                    }
                }
                H[k][j] = H[j][k];
            }
            for (k = 0; k < dF; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[j][dB+k] += dydB[i-j]*dydF[i-k];
                    }
                }
                H[dB+k][j] = H[j][dB+k];
            }
        }

        for (j = 0; j < dF; j++)
        {
            for (k = j; k < dF; k++)
            {
                for (i = 0; i < du; i++)
                {            
                    if ((i>=j)&(i>=k))
                    {
                        H[dB+j][dB+k] += dydF[i-j]*dydF[i-k];
                    }
                }
                H[dB+k][dB+j] = H[dB+j][dB+k];
            }
        }
        

    }
    
    free(dydB);
    free(dydF);
    free(y0);
    free(e);

    free(pB);
    free(pF);

}




void mexFunction(int nlhs, mxArray *plhs[ ],int nrhs, const mxArray *prhs[ ]) 
{
    int i;
    
    du = mxGetM(prhs[0]);
    dteta = mxGetM(prhs[2]);

    u = mxGetPr(prhs[0]);
    ym = mxGetPr(prhs[1]);
    teta = mxGetPr(prhs[2]);
    dB = *mxGetPr(prhs[3]);
    dF = *mxGetPr(prhs[4]);
    delay = *mxGetPr(prhs[5]);


    if (dteta == dB+dF)
    {
        plhs[0] = mxCreateDoubleMatrix(dteta, 1, mxREAL);
        tetafim = mxGetPr(plhs[0]);

        gra = malloc((dteta)*sizeof(double));
        H = malloc(dteta * sizeof(double *));
        H[0] = malloc(dteta * dteta * sizeof(double));
        for(i = 1; i < dteta; i++)
            H[i] = H[0] + i * dteta;

        memcpy(tetafim, teta, sizeof(double) *dteta);
        sib_steepest(tetafim);
        sib_newton(tetafim);
        
        free(gra);
        free(H[0]);
        free(H);

    }else{
        plhs[0] = mxCreateDoubleMatrix(1, 1, mxREAL);
        tetafim = mxGetPr(plhs[0]);
        tetafim = 0;
        
        mexPrintf("ERROR: Length of theta should be nb+nf !!! \n");
    }
    
} 

        

        
        
        
        
