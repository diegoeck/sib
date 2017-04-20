clc
u=(randn(1000,1));
%yi=filter([10 0],[1 -0.9],u);

T=[]
for j=1:10

polos=[0.9 0.8 0.7];

yi=filter(prod(1-polos),poly(polos),u);

%vi=filter([1 -0.9],[1 -0.9],randn(2000,1));
vi=randn(1000,1)*0.01;

y=yi+vi
tao=sib_arx(u,y,3,1,0);
%ta=tao*0.9;
%tic()
%[t2] = sib_oe_conjugate(u,y,ta,0,1); 
%toc()
%tic()
%[t3] = sib_oe_conjugate3(u,y,ta,0,1); 
%toc()

%[teta,J] = sib_oe_hessian(u,y,tao,1000,0,1)
%[teta]=sib_oe(u,y,3,1,0)

%[teta]=sib_oe(u,y,3,1,0)
tarmax=[tao;0;0;0]

%[g,H,J] = sib_armax_grad(u,y,tarmax,0,1,1)

tar = sib_armax_hessian(u,y,tarmax,100,0,1,3)

T=[T tar];
end
%tic()
%ta=sib_arx(u,y,2,2,1);
%ts=sib_sm(u,y,2,2,1);
%to=sib_oe_julia(u,y,ta,100000,1,2);
%to=sib_oe(u,y,2,2,1);

%toc()


% A1=[1;tao(2:4)]';
% B1=[tao(1) 0 0 0];
% C1=[1 0 0 0];
% G1=tf(B1,A1,1);
% H1=tf(C1,A1,1);
% 
% ya=filter([tao(1)],[1;tao(2:4)],u);
% yap=lsim(minreal(inv(H1)*G1),u)+y-lsim(inv(H1),y);
% 
% %ys=filter([0,ts(1),ts(2)],[1;ts(3:4)],u);
% %yo=filter([t3(1)],[1;t3(2:3)],u);
% 
% 
% A2=[1;teta(2:4)]';
% B2=[teta(1) 0 0 0];
% C2=[1 0 0 0];
% G2=tf(B2,A2,1);
% H2=tf(C2,A2,1);
% 
% yh=filter([teta(1)],[1;teta(2:4)],u);
% %yhp=lsim(minreal(inv(H2)*G2),u)+y-lsim(inv(H2),y);
% 
% 
% figure(1)
% plot(y)
% hold on
% plot(ya,'g')
% 
% plot(yh,'r')
% hold off
% 
% figure(2)
% plot(y)
% hold on
% plot(yap,'g')
% 
% plot(yh,'r')
% hold off
% 
