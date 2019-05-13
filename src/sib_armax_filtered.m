function [theta3, m] = sib_armax_filtered(u, y, na, nb, nc, nz);
%  [theta, m] = sib_armax_filtered(u, y, na, nb, nc, nz)
%
%  Prediciton error method with ARMAX structure
%
%         B(z)           C(z)               
%  y(t) = ---- u(t-nz) + ---- e(t) 
%         A(z)           A(z)
%


FF=[.1 .2 .3 .4 .5 .6 .7 .8 .9];

[theta] = sib_arx(u,y,na,nb,nz);
theta = [theta; zeros(nc,1)]

for i=1:length(FF)

    %Filter the data
    [fb,fa] = butter(1,FF(i));
    yf=filter(fb,fa,y);
    uf=filter(fb,fa,u);  

    %Estimate with filtered data

    [theta2] = sib_armax_c(uf,yf,theta,na,nb,nc,nz);
    theta = theta2;

    
end


%Estimate with real data
[theta3] = sib_armax_c(u,y,theta,na,nb,nc,nz);

m.A = [1; theta3(nb+1:nb+na)];
m.B = [zeros(nz,1); theta3(1:nb)];
m.C = [1; theta3(nb+na+1:nb+na+nc)];
m.D = 1;
m.F = 1;


