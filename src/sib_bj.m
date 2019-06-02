function [theta2, m] = sib_bj(u, y, nb, nc, nd, nf, nz)
%  [theta, m] = sib_bj(u, y, nb, nc, nd, nf, nz)
%
%  Prediciton error method with BJ structure
%
%         B(z)           C(z)               
%  y(t) = ---- u(t-nz) + ---- e(t) 
%         F(z)           D(z)
%
% Input : u = input signal
%         y = output signal
%         nb = number of parameter in polynomial B(z) = b_1 + b_2 z^-1 + ... + b_nb z^(1-nb)       
%         nc = number of parameter in polynomial C(z) = 1 + c_1 z^-1 + c_2 z^-2 + ... + c_nc z^-nc       
%         nd = number of parameter in polynomial D(z) = 1 + d_1 z^-1 + d_2 z^-2 + ... + d_nd z^-nd       
%         nf = number of parameter in polynomial F(z) = 1 + f_1 z^-1 + f_2 z^-2 + ... + f_nf z^-nf       
%         nz = input delay
%
% Output : theta = [ b_1 b_2 ... b_nb c_1 c_2 ... c_nc d_1 d_2 ... d_nd f_1 f_2 ... f_nf ]'
%          m = struct with model polynomials 

theta0 = sib_arx(u, y, nf, nb, nz)

if (nd>nf)
    theta1 = [ theta0(nf+1:nb+nf); zeros(nc, 1); theta0(1:nf); zeros(nd-nf, 1); theta0(1:nf) ];
else
    theta1 = [ theta0(nf+1:nb+nf); zeros(nc, 1); theta0(1:nd); theta0(1:nf) ];
end

nb
nc
nd
nf
theta1


theta2 = sib_bj_c(u, y, theta1, nb, nc, nd, nf, nz);

m.A =1;
m.B = [ zeros(nz,1); theta2(1:nb) ];
m.C = [ 1; theta2(nb+1:nb+nc) ];
m.D = [ 1; theta2(nb+nc+1:nb+nc+nd) ];
m.F = [ 1; theta2(nb+nc+nd+1:nb+nc+nd+nf) ];

