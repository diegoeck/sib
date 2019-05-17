function [Torg] = sib_plota(T, a)
% [Torg] = sib_plota(T,a)
%
% Plot Monte-Carlo data T = [ t0 t1 t2 ... tn ], order by parameter *a*.

[b,c] = size(T);
TT = [ T(a,:);T];
TTT = (sortrows([TT']))';
Torg = TTT(2:end,:);
plot(Torg');
xlim([0 c+1]);
xlabel('Monte-Carlo run');
ylabel('\theta');