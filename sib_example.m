
clc

N=1000;
u=randn(N,1);
G=tf([0 1],[1 -0.9],1);
H=tf([1 -0.9],[1 -0.9],1);

yn=lsim(G,u);

T=[];
E=[];

for i=1:1

    y=yn+lsim(H,randn(N,1));

    %[theta] = sib_arx(u,y,1,1,1);
    %[theta]=sib_sm(u,y,1,1,1); 
    %[theta]=sib_oe(u,y,1,1,1); 
    %[theta]=sib_oe_filtered(u,y,1,1,1); 
    [theta]=sib_armax(u,y,1,1,1,1); 
    
    T=[T  theta];

end

sib_plota(T,1)

figure(2)
plot(y)
hold on
plot(yn,'r')
hold off
xlabel('t')
ylabel('y(t)')

