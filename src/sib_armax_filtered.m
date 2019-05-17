function [theta3, m] = sib_armax_filtered(u, y, na, nb, nc, nz);
%  [theta, m] = sib_armax_filtered(u, y, na, nb, nc, nz)
%
%  Prediciton error method with ARMAX structure
%  Improved filtered version with better convergence
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
% Output : theta = [ b_1 b_2 ... b_nb a_1 a_2 ... a_na c_1 c_2 ... c_nc]'
%          m = struct with model polynomials 

FF=[.1 .2 .3 .4 .5 .6 .7 .8 .9];

theta1 = sib_arx(u,y,na,nb,nz);
theta2 = [ theta1; zeros(nc,1) ];

for i=1:length(FF)

    %Filter the data
    [fb,fa] = butter(1, FF(i));
    yf = filter(fb, fa, y);
    uf = filter(fb, fa, u);  

    %Estimate with filtered data

    theta3 = sib_armax_c(uf, yf, theta2, na, nb, nc, nz);
    theta2 = theta3;

end

%Estimate with real data
theta3 = sib_armax_c(u, y, theta2, na, nb, nc, nz);

m.A = [ 1; theta3(nb+1:nb+na) ];
m.B = [ zeros(nz,1); theta3(1:nb) ];
m.C = [ 1; theta3(nb+na+1:nb+na+nc) ];
m.D = 1;
m.F = 1;


