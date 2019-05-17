function [yf] = sib_predict(u, y, m)
% [ys] = sib_predict(u, y, m)
%
% One step ahead prediction
%    
% ys(t) = y(t) + H(z)^(-1) ( G(z)*u(t) - y(t) ) 

yf = filter(m.B,m.A,u);
yf = filter(1,m.F,yf)-y;
yf = filter(m.D,m.C,yf);
yf = y+filter(m.A,1,yf);



