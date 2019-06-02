function [theta2, m] = sib_armax(u, y, na, nb, nc, nz)
%  [theta, m] = sib_armax(u, y, na, nb, nc, nz)
%
%  Prediciton error method with ARMAX structure
%
%         B(z)           C(z)               
%  y(t) = ---- u(t-nz) + ---- e(t) 
%         A(z)           A(z)
%
% Input : u = input signal
%         y = output signal
%         na = number of parameter in polynomial A(z) = 1 + a_1 z^-1 + a_2 z^-2 + ... + a_na z^-na       
%         nb = number of parameter in polynomial B(z) = b_1 + b_2 z^-1 + ... + b_nb z^(1-nb)       
%         nc = number of parameter in polynomial C(z) = 1 + c_1 z^-1 + c_2 z^-2 + ... + c_nc z^-nc       
%         nz = input delay
%
% Output : theta = [ a_1 a_2 ... a_na b_1 b_2 ... b_nb c_1 c_2 ... c_nc]'
%          m = struct with model polynomials 


theta0 = sib_arx(u,y,na,nb,nz);
theta1 = [ theta0; zeros(nc,1) ];
theta2 = sib_armax_c(u, y, theta1, na, nb, nc, nz);

m.A = [ 1; theta2(1:na) ];
m.B = [ zeros(nz,1); theta2(na+1:na+nb) ];
m.C = [ 1; theta2(na+nb+1:na+nb+nc) ];
m.D = 1;
m.F = 1;

