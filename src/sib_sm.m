function [theta2, m] = sib_sm(u, y, nb, nf, nz)
%  [theta, m] = sib_sm(u, y, nb, nf, nz)
%
%  Stieglitz-McBride Method
%
%         B(z)          
%  y(t) = ---- u(t-nz) + e(t) 
%         F(z)        
%
% Input : u = input signal
%         y = output signal
%         nf = number of parameter in polynomial F(z) = 1 + f_1 z^-1 + f_2 z^-2 + ... + f_nf z^-nf       
%         nb = number of parameter in polynomial B(z) = b_1 + b_2 z^-1 + ... + b_nb z^(1-nb)       
%         nz = input delay
%
% Output : theta = [ b_1 b_2 ... b_nb f_1 f_2 ... f_nf ]'
%          m = struct with model polynomials 

theta1 = sib_arx(u, y, nf, nb, nz);

for i=1:100
    uf = filter(1, [ 1; theta1(1:nf) ]', u);
    yf = filter(1, [ 1; theta1(1:nf) ]', y);
    theta1 = sib_arx(uf, yf, nf, nb, nz);
end

theta2 = [ theta1(nf+1:end); theta1(1:nf) ];

m.A = 1;
m.B = [ zeros(nz,1); theta1(nf+1:end) ];
m.C = 1;
m.D = 1;
m.F = [ 1; theta1(1:nf) ];




