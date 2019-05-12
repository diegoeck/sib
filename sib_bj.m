function [theta5]=sib_bj(u,y,na,nb,nc,nd,nz)

[theta3]=sib_armax(u,y,na,nb,nc,nz);
theta3     = [theta3; theta3(nb+1:nb+nc)];

P.u=u;
P.y=y;
P.na=na;
P.nb=nb;
P.nc=nc;
P.nz=nz;

[theta4,J] = sib_gradiente(@sib_bj_grad,theta3,10000,P);
[theta5,J] = sib_newton(@sib_bj_grad,theta4,1000,P);
 
 
