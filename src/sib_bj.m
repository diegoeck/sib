function [theta4, m] = sib_bj(u, y, nf, nb, nc, nd, nz)
%  [theta, m] = sib_bj(u, y, nf, nb, nc, nd, nz)
%
%  Prediciton error method with BJ structure
%
%         B(z)           C(z)               
%  y(t) = ---- u(t-nz) + ---- e(t) 
%         F(z)           D(z)
%
% Input : u = input signal
%         y = output signal
%         nf = number of parameter in polynomial F(z) = 1 + f_1 z^-1 + f_2 z^-2 + ... + f_nf z^-nf       
%         nb = number of parameter in polynomial B(z) = b_1 + b_2 z^-1 + ... + b_nb z^(1-nb)       
%         nc = number of parameter in polynomial C(z) = 1 + c_1 z^-1 + c_2 z^-2 + ... + c_nc z^-nc       
%         nd = number of parameter in polynomial D(z) = 1 + d_1 z^-1 + d_2 z^-2 + ... + d_nd z^-nd       
%         nz = input delay
%
% Output : theta = [ b_1 b_2 ... b_nb f_1 f_2 ... f_nf c_1 c_2 ... c_nc d_1 d_2 ... d_nd ]'
%          m = struct with model polynomials 

theta0 = sib_arx(u, y, nf, nb, nz);
theta1 = [ theta0; zeros(nc,1) ];
%theta2 = sib_armax_c(u, y, theta1, nf, nb, nc, nz);
theta2 = theta1;

if (nd>nf)
    D=[ theta2(nb+1:nb+nf); zeros(nd-nf,1) ];
else
    D=[ theta2(nb+1:nb+nd) ];
end

theta3 = [theta2; D];
theta4 = sib_bj_c(u, y, theta3, nf, nb, nc, nd, nz);

m.A=1;
m.B=[ zeros(nz,1); theta4(1:nb) ];
m.C=[ 1; theta4(nb+nf+1:nb+nf+nc) ];
m.D=[ 1; theta4(nb+nf+nc+1:nb+nf+nc+nd) ];
m.F=[ 1; theta4(nb+1:nb+nf) ];

