function [theta3]=sib_armax(u,y,na,nb,nc,nz)
P.u=u;
P.y=y;
P.na=na;
P.nb=nb;
P.nc=nc;
P.nz=nz;


[theta1]  = sib_arx(u,y,na,nb,nz);

theta1     = [theta1; zeros(nc,1)];

[theta2,J] = sib_gradiente(@sib_armax_grad,theta1,10000,P);
[theta3,J] = sib_newton(@sib_armax_grad,theta2,1000,P);
