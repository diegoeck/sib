
clc

N=1000;
u=randn(N,1);
u=ones(N,1);
G=tf([0 1],[1 -0.9],1);
H=tf([1 0.7],[1 -0.9],1);

yn=lsim(G,u);

T=[];
E=[];

for i=1:1

    y=yn+lsim(H,randn(N,1));

    [theta1, m1] = sib_arx(u, y, 1, 1, 1);
    [theta2, m2] = sib_sm(u, y, 1, 1, 1); 
    [theta3, m3] = sib_oe(u, y, 1, 1, 1); 
    [theta4, m4] = sib_oe_filtered(u, y, 1, 1, 1); 
    [theta5, m5] = sib_armax(u, y, 1, 1, 1, 1); 
    [theta6, m6] = sib_armax_filtered(u, y, 1, 1, 1, 1); 
    
    T=[T  theta1];

end

%sib_plota(T,1)

yf = sib_predict(u, y, m5); 

figure(2)
plot(y)
hold on
plot(yn,'r')
plot(yf,'g')
hold off
xlabel('t')
ylabel('y(t)')

