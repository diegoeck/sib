function [theta2, m] = sib_oe_filtered(u, y, na, nb, nz);
%  [theta, m] = sib_oe_filtered(u, y, nf, nb, nz)
%
%  Prediciton error method with OE structure
%  Improved filtered version with better convergence
%
%         B(z)          
%  y(t) = ---- u(t-k) + e(t) 
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

FF=[.1 .2 .3 .4 .5 .6 .7 .8 .9];

theta1 = sib_arx(u,y,na,nb,nz);

for i = 1:length(FF)

    %Filter the data
    [fb,fa] = butter(1,FF(i));
    yf = filter(fb, fa, y);
    uf = filter(fb, fa, u);  

    %Estimate with filtered data
    theta2 = sib_oe_c(uf, yf, theta1, nb, nz) ;
    theta1 = theta2;

end


%Estimate with real data
theta2 = sib_oe_c(u,y,theta,nb,nz) ;

m.A = 1;
m.B = [ zeros(nz,1); theta2(1:nb) ];
m.C = 1;
m.D = 1;
m.F = [ 1; theta2(nb+1:end) ];


