#include <stdint.h>
#include <math.h>
#include <string.h>
#include <float.h>
#include "mex.h"
#include "lapack.h"
#include <sib_basic.h>
#include <sib_oe_c.h>
typedef uint16_t char16_t;

extern double *teta;
extern double *gra;
extern double **H;
extern int dteta;
extern int mode;

extern bool utIsInterruptPending();


void sib_steepest(double* tetafim)
{

    int i,j,k;
    double *tetatemp;
    double *dir;
    double passo = 0.00001;
    double N;
    double J1[1];
    double J2[1];

    tetatemp=malloc((dteta)*sizeof(double));
    dir=malloc((dteta)*sizeof(double));
    
    memcpy(tetafim,teta,sizeof(double) *dteta);

    for(i=0; i<100; i++)
    {

        for(j=0; j<100; j++)
        {   
            mode=1;
            grad(tetafim,J1);
    
            if((i==0)&(j==0))
            {
                memcpy(dir,gra,sizeof(double)*dteta);
            }

            for (k = 0; k < dteta; k++)
            {
                dir[k] = (4*dir[k]+gra[k])/5;
            }
    
            N = norm(dir,dteta); 

            for (k = 0; k < dteta; k++)
            {
                tetatemp[k] = tetafim[k] - passo*dir[k]/N;
            }

            mode=0;
            grad(tetatemp,J2);

            if (J2[0]>J1[0])
            {
                passo=passo*0.99;

            }else{   
                passo=passo*1.01;
                memcpy(tetafim,tetatemp,sizeof(double) *dteta);
            }
    
        }    
    
        
        mexPrintf("G: ");
        for (k = 0; k < (int)floor(i/5); k++)
        {
            mexPrintf("#");
        }
        for (k = (int)floor(i/5); k < 20; k++)
        {
            mexPrintf("-");
        }
        mexPrintf(" %2d%% %1.10f %1.10f \n",i,J1[0],passo);
        
        if (utIsInterruptPending()) {       
            mexPrintf("Ctrl-C Detected. END\n\n");
            return;
        }
        
        mexEvalString("drawnow;");

        

        if (passo<0.0000001)
        {
            break;
        }    
    }

    free(tetatemp);
    free(dir);

}
 



void sib_newton(double* tetafim)
{

    int i,j,k;
    double *tetatemp;
    double passo;
    double J1[1];

    long Ns = (long)dteta;
    long Ms = 1;

    long *ipiv;
    long info;

    ipiv=malloc((dteta)*sizeof(long));

    tetatemp=malloc((dteta)*sizeof(double));

    memcpy(tetafim,teta,sizeof(double) *dteta);


    for(i=0;i<100;i++)
    {
        mode=2;
        grad(tetafim,J1);

        //mexPrintf("G: \n");
        //mexPrintf("%1.10f %1.10f %1.10f\n",gra[0],gra[1],gra[2]);
        
        //mexPrintf("H: \n");
        //mexPrintf("%1.10f %1.10f %1.10f %1.10f\n",H[0][0],H[0][1],H[0][2],H[0][3]);
        //mexPrintf("%1.10f %1.10f %1.10f %1.10f\n",H[1][0],H[1][1],H[1][2],H[1][3]);
        //mexPrintf("%1.10f %1.10f %1.10f %1.10f\n",H[2][0],H[2][1],H[2][2],H[2][3]);
        //mexPrintf("%1.10f %1.10f %1.10f %1.10f\n",H[3][0],H[3][1],H[3][2],H[3][3]);        

        //mexPrintf("%1.10f %1.10f %1.10f\n",H[0][0],H[0][1],H[0][2]);
        //mexPrintf("%1.10f %1.10f %1.10f\n",H[1][0],H[1][1],H[1][2]);
        //mexPrintf("%1.10f %1.10f %1.10f\n",H[2][0],H[2][1],H[2][2]);

        
        dgesv(&Ns, &Ms, *H, &Ns, ipiv, gra, &Ns, &info);
        //mexPrintf("%1.10f %1.10f %1.10f\n",gra[0],gra[1],gra[2]);

    
        for (k = 0; k < dteta; k++)
        {
            tetatemp[k] = tetafim[k] - gra[k]*i/100;
        }
    
        memcpy(tetafim,tetatemp,sizeof(double) *dteta);
    
        
        mexPrintf("N: ");
        for (k = 0; k < (int)floor(i/5); k++)
        {
            mexPrintf("#");
        }
        for (k = (int)floor(i/5); k < 20; k++)
        {
            mexPrintf("-");
        }
        mexPrintf(" %02d%% %1.10f \n",i,J1[0]);
        
        
        
        if (utIsInterruptPending()) {       
            mexPrintf("Ctrl-C Detected. END\n\n");
            return;
        }
        
        mexEvalString("drawnow;");
        
    }

    free(tetatemp);

}
 


