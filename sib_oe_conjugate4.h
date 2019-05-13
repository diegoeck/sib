#ifndef OE_H
#define OE_H

void cost(double* teta, double* J);
void grad(double* teta, int dteta, int delay, int dnum, double* u, int du,double* ym, double* J,double* gra);

#endif
