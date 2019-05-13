function [theta2, m] = sib_oe(u, y, nf, nb, nz)
%  [theta, m] = sib_oe(u, y, nf, nb, nz)
%
%  Prediciton error method with OE structure
%
%         B(z)          
%  y(t) = ---- u(t-nz) + e(t) 
%         F(z)        
%

[theta] = sib_arx(u,y,nf,nb,nz);
[theta2] = sib_oe_c(u,y,theta,nb,nz);


m.A=1;
m.B=[zeros(nz,1); theta2(1:nb)];
m.C=1;
m.D=1;
m.F=[1; theta2(nb+1:end)];




