FF=[.9 .8 .7 .6 .5 .4 .3 .2 .1]*40;
FF=exp(log(0.05)./FF)


figure(1)
hold on
figure(2)
hold on
figure(3)
hold on

for i=1:9

G=zpk([],[FF(i)],1,1)
G=G/dcgain(G);

figure(1)
step(G)

figure(2)
impulse(G)

figure(3)
bode(G)
end
