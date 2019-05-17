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
int dD;
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
    double *y4;
    double *y5;
    double *y6;
    double *y7;

    int i,j,k;

    y0=malloc((du)*sizeof(double));
    y1=malloc((du)*sizeof(double));
    y2=malloc((du)*sizeof(double));
    y3=malloc((du)*sizeof(double));
    y4=malloc((du)*sizeof(double));
    y5=malloc((du)*sizeof(double));
    y6=malloc((du)*sizeof(double));
    y7=malloc((du)*sizeof(double));
    e=malloc((du)*sizeof(double));

    
    //mexPrintf(" %f %f %f %f \n",teta[0],teta[1],teta[2],teta[2]);
    
    A=malloc((dB+dA)*sizeof(double));
    for (k = 0; k < dB+dA; k++)
    {
        A[k]=teta[k];
    }
    
    filtra(A, dB, dA, delay, u, du, y0); 
    free(A);
    
    for (k = 0; k < du; k++)
    {
        y0[k]-=ym[k];
    }
    
    A=malloc((dC+dD+1)*sizeof(double));
    A[0]=1;
    for (k = 0; k < dD; k++)
    {
        A[k+1]=teta[dB+dA+dC+k];
    }
    for (k = 0; k < dC; k++)
    {
        A[k+1+dD]=teta[dB+dA+k];
    }

    filtra(A, dD+1, dC, 0, y0, du, y1);    
    free(A);
    
    J[0]=0;
    for (k = 0; k < du; k++)
    {
        e[k]=(y1[k]);
        J[0]+=(e[k])*(e[k]);
    }

    
    
    if(mode==1 | mode==2){
        
          
        for (j = 0; j < dA+dB+dC+dD; j++)
        {
            gra[j]=0;
        }
      
        // Calcula gradiente com relacao a B
        A=malloc((dC+dD+1)*sizeof(double));
        A[0]=1;
        for (k = 0; k < dD; k++)
        {
            A[k+1]=teta[dB+dA+dC+k];
        }
        for (k = 0; k < dC; k++)
        {
            A[k+1+dD]=teta[dB+dA+k];
        }
        filtra(A, dD+1, dC, 0, u, du, y1);    
        free(A);

        A=malloc((1+dA)*sizeof(double));
        A[0]=1;
        for (k = 0; k < dA; k++)
        {
            A[k+1]=teta[dB+k];
        }
        filtra(A, 1, dA, delay, y1, du, y2);
        free(A);
        

        for (j = 0; j < dB; j++)
        {
            gra[j]=0;
            for (k = j; k < du; k++)
            {
                gra[j]+=(e[k])*(y2[k-j]);
            }
        }

        // Calcula gradiente com relacao a A
        A=malloc((dB+dA)*sizeof(double));
        for (k = 0; k < dB+dA; k++)
        {
            A[k]=teta[k];
        }
        
        filtra(A, dB, dA, delay, y1, du, y3);    
        free(A);

        A=malloc((1+dA)*sizeof(double));
        A[0]=-1;
        for (k = 0; k < dA; k++)
        {
            A[k+1]=teta[dB+k];
        }
        filtra(A, 1, dA, 1, y3, du, y4);
        free(A);
        

        for (j = 0; j < dA; j++)
        {
            gra[dB+j]=0;
            for (k = j; k < du; k++)
            {
                gra[dB+j]+=(e[k])*(y4[k-j]);
            }
        }
        
        
        // Calcula gradiente com relacao a C
        A=malloc((dC+dD+1)*sizeof(double));
        A[0]=1;
        for (k = 0; k < dD; k++)
        {
            A[k+1]=teta[dB+dA+dC+k];
        }
        for (k = 0; k < dC; k++)
        {
            A[k+1+dD]=teta[dB+dA+k];
        }
        filtra(A, dD+1, dC, 0, y0, du, y5);    
        free(A);
        
        A=malloc((dC+1)*sizeof(double));
        A[0]=-1;
        for (k = 0; k < dC; k++)
        {
            A[k+1]=teta[dB+dA+k];
        }
        filtra(A, 1, dC, 1, y5, du, y6);    
        free(A);
        
        for (j = 0; j < dC; j++)
        {
            gra[dA+dB+j]=0;
            for (k = j; k < du; k++)
            {
                gra[dA+dB+j]+=(e[k])*(y6[k-j]);
            }
        }
        
        // Calcula gradiente com relacao a D
        A=malloc((dC+1)*sizeof(double));
        A[0]=1;
        for (k = 0; k < dC; k++)
        {
            A[k+1]=teta[dB+dA+k];
        }
        filtra(A, 1, dC, 1, y0, du, y7);    
        free(A);
        
        for (j = 0; j < dD; j++)
        {
            gra[dA+dB+dC+j]=0;
            for (k = j; k < du; k++)
            {
                gra[dA+dB+dC+j]+=(e[k])*(y7[k-j]);
            }
        }

    }
    
    //mexPrintf("Teta %1.10f %1.10f %1.10f %1.10f\n",gra[0],gra[1],gra[2],gra[3]);

    
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
                        H[j][k]+=(y2[i-k])*(y2[i-j]);
                        H[k][j]=H[j][k];
                    }
                }
                for (k = dB; k < dB+dA; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                        H[j][k]+=(y2[i-k])*(y4[i-j]);
                        H[k][j]=H[j][k];
                    }
                }

                for (k = dB+dA; k < dB+dA+dC; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                        H[j][k]+=(y2[i-j])*(y6[i-k]);
                        H[k][j]=H[j][k];
                    }
                }
                for (k = dB+dA+dC; k < dB+dA+dC+dD; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                        H[j][k]+=(y2[i-j])*(y7[i-k]);
                        H[k][j]=H[j][k];
                    }
                }

            }
        

            for (j = dB; j < dB+dA; j++)
            {
                for (k = j; k < dB+dA; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                        H[j][k]+=(y4[i-k])*(y4[i-j]);
                        H[k][j]=H[j][k];
                    }
                }
                for (k = dB+dA; k < dB+dA+dC; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                        H[j][k]+=(y4[i-j])*(y6[i-k]);
                        H[k][j]=H[j][k];
                    }
                }
                for (k = dB+dA+dC; k < dB+dA+dC+dD; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                        H[j][k]+=(y4[i-j])*(y7[i-k]);
                        H[k][j]=H[j][k];
                    }
                }

            }
            
            for (j = dB+dA; j < dB+dA+dC; j++)
            {
                for (k = j; k < dB+dA+dC; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                     H[j][k]+=(y6[i-j])*(y6[i-k]);
                     H[k][j]=H[j][k];
                    }
                }
                for (k = dB+dA+dC; k < dB+dA+dC+dD; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                     H[j][k]+=(y6[i-j])*(y7[i-k]);
                     H[k][j]=H[j][k];
                    }
                }

            }
    
            
            for (j = dB+dA+dC; j < dB+dA+dC+dD; j++)
            {
                for (k = j; k < dB+dA+dC+dD; k++)
                {
                    if( (i>=j) & (i>=k) ) 
                    {
                     H[j][k]+=(y7[i-j])*(y7[i-k]);
                     H[k][j]=H[j][k];
                    }
                }

            }
    
        
        }

    }
    
    free(y);
    free(y0);
    free(y1);
    free(y2);
    free(y3);
    free(y4);
    free(y5);
    free(y6);
    free(y7);
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
    dD = *mxGetPr(prhs[6]);
    delay = *mxGetPr(prhs[7]);
    
    

    gra=malloc((dteta)*sizeof(double));

    H = malloc(dteta * sizeof(double *));
    H[0] = malloc(dteta * dteta * sizeof(double));
    for(i = 1; i < dteta; i++)
        H[i] = H[0] + i * dteta;


    plhs[0] = mxCreateDoubleMatrix(dteta,1,mxREAL);
    tetafim = mxGetPr(plhs[0]);

    //mexPrintf(" %d %d %d %d  \n",dteta,dA,dB,dC);

    //mexPrintf("Teta %1.10f %1.10f %1.10f %1.10f\n",teta[0],teta[1],teta[2],teta[3]);

    sib_steepest(tetafim);
    memcpy(teta,tetafim,sizeof(double) *dteta);
    sib_newton(tetafim);
    
    free(gra);
    free(H);
} 

        

        
        
        
        
