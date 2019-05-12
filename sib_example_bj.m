clc
close all

N=1000;
u=randn(N,1);
u=ones(N,1);
u(1:100)=zeros(100,1);
G=tf(zpk([ ],[0.9],1,1));
H=tf(zpk([0.7],[.6],1,1));

u=randn(N,1);
yn=lsim(G,u);
y=yn+lsim(H,randn(N,1))*0.1;

theta_armax=sib_bj(u,y,1,1,1,1,1); 

%Z=iddata(y,u)


GG=tf([0 theta_armax(1)],[1 theta_armax(2)],1) 
HH=tf([1 theta_armax(3)],[1 theta_armax(4)],1) 

%GG=tf([0 0 theta_armax(1)],[1 theta_armax(2) theta_armax(3)],1) 
%HH=tf([1 theta_armax(4) theta_armax(5)],[1 theta_armax(6) theta_armax(7)],1) 


figure
plot(y)
hold on
yp1=lsim(inv(HH),lsim(GG,u)-y)+y;
yp2=lsim(inv(H),lsim(G,u)-y)+y;
plot(yp1)
plot(yp2)
hold off