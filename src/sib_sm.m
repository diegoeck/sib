function [theta, m] = sib_sm(u, y, nf, nb, nz)
%  [theta, m] = sib_sm(u, y, nf, nb, nz)
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

theta = sib_arx(u, y, nf, nb, nz);

for i=1:10
    uf = filter(1, [ 1; theta(nb+1:nb+nf) ]', u);
    yf = filter(1, [ 1; theta(nb+1:nb+nf) ]', y);
    theta = sib_arx(uf, yf, nf, nb, nz);
end

m.A = 1;
m.B = [ zeros(nz,1); theta(1:nb) ];
m.C = 1;
m.D = 1;
m.F = [ 1; theta(nb+1:end) ];




