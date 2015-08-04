
clc

N=500;
u=randn(N,1);
G=tf(zpk([0.9],[.85 .95],1,1));
yn=lsim(G,u);

T=[];
E=[];

for i=1:100
i
    y=yn+randn(N,1)*1;

%      [theta]=sib_oe(u,y,2,2,1); %98 %60

% %%% OE 71     48
%   dat=iddata(y,u,1);
%   Model = oe(dat,[2 2 1]);
%   [mA,mB,mC,mD,mF] = polydata(Model);
%   theta=[mB(2:3) mF(2:3)]';
% % 

     [theta,estavel]=sib_oe_filtered(u,y,2,2,1); %100% 2

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

