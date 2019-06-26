clc
clear all

G = tf([zeros(1,20) .1], [1 -0.9 zeros(1,20)], 1);
H = tf([1 0 0], [1 -1.9 0.98], 1);
H = tf([1 0 0], [1 0 0], 1);

C=1/(1-G)

%C = tf([1 -0.8]/3, [1 -1], 1);



T = feedback(C*G,1)
S = feedback(1,C*G)


N = 1000;
r = square((1:N)*(2*pi)/1000)';
%v = lsim(H,randn(N,1))/100;
v=sin(0.075*(1:N))/10;

y = lsim(T,r) + lsim(S,v);
u = lsim(C*S,r) - lsim(C*S,v);

[theta, m] = sib_bj(u, y, 1, 0, 2, 1, 21)

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
