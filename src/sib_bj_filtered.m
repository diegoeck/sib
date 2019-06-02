function [theta3, m] = sib_bj_filtered(u, y, nb, nc, nd, nf, nz);
%  [theta, m] = sib_bj_filtered(u, y, nb, nc, nd, nf, nz)
%
%  Prediciton error method with BJ structure
%  Improved filtered version with better convergence
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

FF=[.1 .2 .3 .4 .5 .6 .7 .8 .9];

theta0 = sib_arx(u, y, nf, nb, nz);
if (nd>nf)
    theta1 = [ theta0(nb+1:nb+nf); zeros(nc, 1); theta0(1:nf); zeros(nd-nf, 1); theta0(1:nf) ];
else
    theta1 = [ theta0(nb+1:nb+nf); zeros(nc, 1); theta0(1:nf); theta0(1:nd) ];
end

for i=1:length(FF)

    %Filter the data
    [fb,fa] = butter(1, FF(i));
    yf = filter(fb, fa, y);
    uf = filter(fb, fa, u);  

    %Estimate with filtered data
    theta2 = sib_bj_c(uf, yf, theta1, nb, nc, nd, nf, nz);
    theta1 = theta2;

    
end

%Estimate with real data
theta3 = sib_bj_c(u, y, theta1, nb, nc, nd, nf, nz);

m.A = 1;
m.B = [ zeros(nz,1); theta3(1:nb) ];
m.C = [ 1; theta3(nb+1:nb+nc) ];
m.D = [ 1; theta3(nb+nc+1:nb+nc+nd) ];
m.F = [ 1; theta3(nb+nc+nd+1:nb+nc+nd+nf) ];
