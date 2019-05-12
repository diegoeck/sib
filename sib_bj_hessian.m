function [teta,J] = sib_bj_hessian(u,y,teta,NN,nz,na,nb,nc)

P.u=u;
P.y=y;
P.na=na;
P.nb=nb;
P.nc=nc;
P.nz=nz;

[teta,J]=sib_gradiente(teta,1000,P)
[teta,J]=sib_newton(teta,100,P)
