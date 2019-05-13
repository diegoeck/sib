function [theta, m] = sib_sm(u, y, na, nb, nz)
%  [theta, m] = sib_sm(u, y, na, nb, nz)
%
%  Stieglitz-McBride Method
%
%         B(z)          
%  y(t) = ---- u(t-nz) + e(t) 
%         F(z)        
%

theta = sib_arx(u,y,na,nb,nz);

for i=1:10
    uf = filter(1, [1; theta(nb+1:nb+na)]', u);
    yf = filter(1, [1; theta(nb+1:nb+na)]', y);
    theta = sib_arx(uf, yf, na, nb, nz);
end

m.A=1;
m.B=[zeros(nz,1); theta(1:nb)];
m.C=1;
m.D=1;
m.F=[1; theta(nb+1:end)];




