function [theta2]=sib_armax(u,y,na,nb,nc,nz)

[theta0] = sib_arx(u,y,na,nb,nz);
[theta1] = [theta0; zeros(nc,1)]
[theta2] = sib_armax_c(u,y,theta1,na,nb,nc,nz)

