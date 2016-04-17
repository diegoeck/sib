
clc

N=1000;
u=randn(N,1);
%u=[0; ones(N-1,1)];
G=tf(zpk([0.85 .75 ],[.7 .8 .9],1,1));


%G=tf([0 1 0.5 -2 1],[1 -1.5 0.7 0.3 -0.2],1);
%H=tf([1 -0.6 0.4],[1 -1.95 0.9506],1);






yn=lsim(G,u);

T=[];

for j=1:100

%v=randn(N,1)*sqrt(norm(G)^2/100);

v=randn(N,1);
%v=lsim(H,v);
v=v*sqrt(var(yn)/var(v)/10000*5);
%v=v*sqrt(var(yn)/var(v)/10000);

y=yn+v;

nb=45;
na=45;
[thetai] = sib_arx(u,y,na,nb,1);
F=tf([1; thetai(nb+1:nb+na)]',[1 zeros(1,na)],1);
uh=lsim(F,u);
yh=lsim(F,y);

max(abs(roots([1; thetai(nb+1:nb+na)]')))


na=3;
nb=3;
[theta]=sib_oe_filtered(u,y,na,nb,1) %98 %60

T=[T theta];

end


figure(1)
plot(y)
hold on
plot(yn,'r')
hold off

abs(roots([1; theta(na+1:na+nb)]'))