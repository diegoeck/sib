clc
clear all

N = 100;
u = square((1:N)*(2*pi)/100)';

G = zpk([], [0.9 0.8 0.7], 1, 1);
G=G/dcgain(G);
H = tf([1], [1], 1);


P1=[]
P2=[]
P3=[]
for i=1:100

    y = lsim(G,u) + lsim(H,randn(N,1))/10;

    Z=iddata(y,u);
    M=oe(Z,[1 3 3])
    P1=[P1 [M.B';M.F']];

    [theta, m] = sib_oe(u, y, 1, 3, 3);
    P2=[P2 theta];

    [theta, m] = sib_oe_filtered(u, y, 1, 3, 3);
    P3=[P3 theta];

    
end

yp = sib_predict(u, y, m); 
ys = sib_simulate(u, m); 
 
figure(1)
plot(y, 'bo-')
hold on
plot(yp, 'go-')
plot(ys, 'ro-')
hold off
xlabel('t')
ylabel('y(t)')
legend('Data', 'Prediction', 'Simulation') 

figure(2)
plot(P1)

figure(3)
plot(P2)

figure(4)
plot(P3)
