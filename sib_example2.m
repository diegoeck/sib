
clc

N=500;
u=randn(N,1);
G=tf(zpk([0.9],[.85 .95],1,1));
yn=lsim(G,u);

H=tf(zpk([0.8],[0.7],1,1));


T=[];
E=[];




for i=1:100
i
erro=lsim(H,randn(N,1)*1);
y=yn+erro;

      [theta]=sib_oe(u,y,2,2,1); %98 %60

 %%% OE 71     48
%   dat=iddata(y,u,1);
%   Model = oe(dat,[2 2 1]);
%   [mA,mB,mC,mD,mF] = polydata(Model);
%   theta=[mB(2:3) mF(2:3)]';
% % 

%     [theta,estavel]=sib_oe_filtered(u,y,2,2,1); %100% 2

    T=[T  theta];

    figure(1)
    sib_plota(T,3)

end

 
figure(2)
plot(y)
hold on
plot(yn,'r')
hold off
xlabel('t')
ylabel('y(t)')

BINS=50
figure(3)
subplot(2,2,1)
hist(T(1,:),BINS)
subplot(2,2,2)
hist(T(2,:),BINS)
subplot(2,2,3)
hist(T(3,:),BINS)
subplot(2,2,4)
hist(T(4,:),BINS)


