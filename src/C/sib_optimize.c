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
    
    //memcpy(tetafim,teta,sizeof(double) *dteta);

    for(i=1; i<=100; i++)
    {

        

        for(j=1; j<=100; j++)
        {   
            mode=1;
            grad(tetafim,J1);
            
            //mexPrintf("T %1.10f %1.10f %1.10f %1.10f\n",tetafim[0],tetafim[1],tetafim[2],tetafim[3]);
            //mexPrintf("G %1.10f %1.10f %1.10f %1.10f\n",gra[0],gra[1],gra[2],gra[3]);

            if((i==1)&(j==1))
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
        for (k = 0; k < i/5; k++)
        {
            mexPrintf("#");
        }
        for (k = i/5; k < 20; k++)
        {
            mexPrintf("-");
        }
        mexPrintf(" %03d%% %1.10f %1.10f \n", i, J1[0], passo);

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
    double passo = 0.00001;
    double J1[1];
    double J2[1];

    //#ifdef _WIN64
    //long long Ns = (long long)dteta;
    //long long Ms = 1;
    //long long *ipiv;
    //ipiv=malloc((dteta)*sizeof(long long));
    //long long info;
    //#else
    ptrdiff_t Ns = dteta;
    ptrdiff_t Ms = 1;
    ptrdiff_t info = 0;
    
    char uplo[] = "U";
    //ptrdiff_t *ipiv;
    //ipiv=malloc((dteta)*sizeof(ptrdiff_t));

    //#endif
    

    
    //int erro=0;


    tetatemp=malloc((dteta)*sizeof(double));

    //memcpy(tetafim,teta,sizeof(double) *dteta);


    for(i=1; i<=1000; i++)
    {

        mode=2;
        grad(tetafim,J1);

        //mexPrintf("G: \n");
        //mexPrintf("antes %1.10f %1.10f %1.10f %1.10f\n",gra[0],gra[1],gra[2],gra[3]);
        
        //mexPrintf("antes: %f\n", gra[3]);
        //mexPrintf("H: \n");
        //mexPrintf("%1.10f %1.10f %1.10f %1.10f\n",H[0][0],H[0][1],H[0][2],H[0][3]);
        //mexPrintf("%1.10f %1.10f %1.10f %1.10f\n",H[1][0],H[1][1],H[1][2],H[1][3]);
        //mexPrintf("%1.10f %1.10f %1.10f %1.10f\n",H[2][0],H[2][1],H[2][2],H[2][3]);
        //mexPrintf("%1.10f %1.10f %1.10f %1.10f\n",H[3][0],H[3][1],H[3][2],H[3][3]);        

        //mexPrintf("%1.10f %1.10f %1.10f\n",H[0][0],H[0][1],H[0][2]);
        //mexPrintf("%1.10f %1.10f %1.10f\n",H[1][0],H[1][1],H[1][2]);
        //mexPrintf("%1.10f %1.10f %1.10f\n",H[2][0],H[2][1],H[2][2]);

        
        
        //dgesv(&Ns, &Ms, *H, &Ns, ipiv, gra, &Ns, &info);
        dposv(uplo, &Ns, &Ms, *H, &Ns, gra, &Ns, &info);

        //mexPrintf("depois %1.10f %1.10f %1.10f %1.10f\n",gra[0],gra[1],gra[2],gra[3]);
        
        //mexPrintf("info: %d\n", info);

    /*
        if (isnan(gra[0]))
        {
            break;
        }
     */
        
        for (k = 0; k < dteta; k++)
        {
            tetatemp[k] = tetafim[k] - gra[k]*(i+1)/1000;
        }
    
        //memcpy(tetafim,tetatemp,sizeof(double) *dteta);
        
        
        for (j = 0; j < 10; j++)
        {
    
            mode=0;
            grad(tetatemp,J2);

            if (J2[0]>J1[0] || ~isfinite(tetatemp[0]))
            {
                //mexPrintf("ERRO \n");
                for (k = 0; k < dteta; k++)
                {
                    tetatemp[k] = (tetatemp[k]+tetafim[k])/2;
                }

            }else{   
                
                memcpy(tetafim,tetatemp,sizeof(double) *dteta);
                break;
                
            }
        }

        if (j==10)
        {
            break;
        }

        
        
        
        if (i%10==0)
        {
            mexPrintf("N: ");
            for (k = 0; k < (i/50); k++)
            {
                mexPrintf("#");
            }
            for (k = (i/50); k < 20; k++)
            {
                mexPrintf("-");
            }
            mexPrintf(" %03d%% %1.10f %1.10f \n", i/10, J1[0], ((float)i)/1000.0);
        }
        
        
        
        if (utIsInterruptPending()) {       
            mexPrintf("Ctrl-C Detected. END\n\n");
            return;
        }
        
        mexEvalString("drawnow;");
        
    }

    //free(ipiv);
    free(tetatemp);

}
 


