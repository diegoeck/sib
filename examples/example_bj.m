clc
clear all

N = 1000;
u = square((1:N)*(2*pi)/100)';
G = tf([0 1],[1 -0.9],1);
H=tf([1 -0.5 0],[1 -1.9 0.98],1);

y = lsim(G,u) + lsim(H,randn(N,1))/10;

[theta, m] = sib_bj(u, y, 1, 2, 2, 2, 1)

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
