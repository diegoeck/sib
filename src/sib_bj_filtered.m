function [theta4, m] = sib_bj_filtered(u, y, na, nb, nc, nd, nz);
%  [theta, m] = sib_bj_filtered(u, y, na, nb, nc, nd, nz)
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
%         nf = number of parameter in polynomial F(z) = 1 + f_1 z^-1 + f_2 z^-2 + ... + f_nf z^-nf       
%         nb = number of parameter in polynomial B(z) = b_1 + b_2 z^-1 + ... + b_nb z^(1-nb)       
%         nc = number of parameter in polynomial C(z) = 1 + c_1 z^-1 + c_2 z^-2 + ... + c_nc z^-nc       
%         nd = number of parameter in polynomial D(z) = 1 + d_1 z^-1 + d_2 z^-2 + ... + d_nd z^-nd       
%         nz = input delay
%
% Output : theta = [ b_1 b_2 ... b_nb f_1 f_2 ... f_nf c_1 c_2 ... c_nc d_1 d_2 ... d_nd ]'
%          m = struct with model polynomials 

FF=[.1 .2 .3 .4 .5 .6 .7 .8 .9];

theta0 = sib_arx(u, y, na, nb, nz);
theta1 = [ theta0; zeros(nc,1) ];
%theta2 = sib_armax_c(u, y, theta1, na, nb, nc, nz);
theta2 = theta1;
if (nd>na)
    D=[ theta2(nb+1:nb+na); zeros(nd-na,1) ];
else
    D=[ theta2(nb+1:nb+nd) ];
end
theta3 = [ theta2; D ];

for i=1:length(FF)

    %Filter the data
    [fb,fa] = butter(1, FF(i));
    yf = filter(fb, fa, y);
    uf = filter(fb, fa, u);  

    %Estimate with filtered data
    theta4 = sib_bj_c(uf, yf, theta3, na, nb, nc, nd, nz);
    theta3 = theta4;

    
end

%Estimate with real data
theta4 = sib_bj_c(u, y, theta3, na, nb, nc, nd, nz);

m.A =1;
m.B = [ zeros(nz,1); theta4(1:nb) ];
m.C = [ 1; theta4(nb+na+1:nb+na+nc) ];
m.D = [ 1; theta4(nb+na+nc+1:nb+na+nc+nd) ];
m.F = [ 1; theta4(nb+1:nb+na) ];



