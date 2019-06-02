function [theta2, m] = sib_oe(u, y, nb, nf, nz)
%  [theta, m] = sib_oe(u, y, nf, nb, nz)
%
%  Prediciton error method with OE structure
%
%         B(z)          
%  y(t) = ---- u(t-nz) + e(t) 
%         F(z)        
%
% Input : u = input signal
%         y = output signal
%         nb = number of parameter in polynomial B(z) = b_1 + b_2 z^-1 + ... + b_nb z^(1-nb)       
%         nf = number of parameter in polynomial F(z) = 1 + f_1 z^-1 + f_2 z^-2 + ... + f_nf z^-nf       
%         nz = input delay
%
% Output : theta = [ b_1 b_2 ... b_nb f_1 f_2 ... f_nf ]'
%          m = struct with model polynomials 

theta0 = sib_arx(u, y, nf, nb, nz)
theta1 = [ theta0(nf+1:nb+nf); theta0(1:nf) ]
theta2 = sib_oe_c(u, y, theta1, nb, nf, nz)

m.A = 1;
m.B = [ zeros(nz,1); theta2(1:nb) ];
m.C = 1;
m.D = 1;
m.F = [ 1; theta2(nb+1:end) ];




