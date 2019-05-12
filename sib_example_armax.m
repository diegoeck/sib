
clc

N=10000;
u=randn(N,1);
u=ones(N,1);
u(1:100)=zeros(100,1)
G=tf(zpk([0 ],[.8 .9],1,1));
H=tf(zpk([0.7 0.6],[.8 .9],1,1));


yn=lsim(G,u);

T=[];
E=[];


y=yn+lsim(H,randn(N,1))*1;
%y=yn+randn(N,1)*1;

theta_armax=sib_armax(u,y,2,1,2,1); 
%theta_oe=sib_oe(u,y,2,1,1)
%theta=sib_armax_filtered(u,y,2,2,2,1); 

Z=iddata(y,u)
M=armax(Z,[2 1 2 1])
%M2=oe(Z,[1 2 1])


