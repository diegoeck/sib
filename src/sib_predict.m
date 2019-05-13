function [yf] = sib_predict(u, y, m)

yf=filter(m.B,m.A,u);
yf=filter(1,m.F,yf)-y;
yf=filter(m.D,m.C,yf);
yf=y+filter(m.A,1,yf);



