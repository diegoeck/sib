function [theta2, m] = sib_armax(u, y, na, nb, nc, nz)
%  [theta, m] = sib_armax(u, y, na, nb, nc, nz)
%
%  Prediciton error method with ARMAX structure
%
%         B(z)           C(z)               
%  y(t) = ---- u(t-nz) + ---- e(t) 
%         A(z)           A(z)
%

[theta0] = sib_arx(u,y,na,nb,nz);
[theta1] = [theta0; zeros(nc,1)]
[theta2] = sib_armax_c(u,y,theta1,na,nb,nc,nz)

m.A=[1; theta2(nb+1:nb+na)];
m.B=[zeros(nz,1); theta2(1:nb)];
m.C=[1; theta2(nb+na+1:nb+na+nc)];
m.D=1;
m.F=1;

