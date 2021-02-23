function [theta, m] = sib_arx(u, y, na, nb, nz)
%  [theta, m] = sib_arx(u, y, na, nb, nz)
%
%  Prediciton error method with ARX structure
%
%         B(z)            1
%  y(t) = ---- u(t-nz) + ---- e(t) 
%         A(z)           A(z)
%
% Input : u = input signal
%         y = output signal
%         na = number of parameter in polynomial A(z) = 1 + a_1 z^-1 + a_2 z^-2 + ... + a_na z^-na       
%         nb = number of parameter in polynomial B(z) = b_1 + b_2 z^-1 + ... + b_nb z^(1-nb)       
%         nz = input delay
%
% Output : theta = [ a_1 a_2 ... a_na  b_1 b_2 ... b_nb ]'
%          m = struct with model polynomials 

phi = [];

for i = 1:na
    phi = [ phi [zeros(i,1) ; -y(1:end-i)] ];
end

for i = 1+nz:nb+nz
    phi = [ phi [zeros(i-1,1); u(1:end-i+1)] ];
    
end

retirar= 1 + max(na,nb+nz-1); %Remove ZEROS iniciais

theta = phi(retirar:end,:)\y(retirar:end);

m.A = [ 1; theta(1:na) ];
m.B = [ zeros(nz,1); theta(na+1:end) ];
m.C = 1;
m.D = 1;
m.F = 1;
