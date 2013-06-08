function [Torg]=sib_plota(T,a)

[b,c]=size(T);
TT=[ T(a,:);T];
TTT=(sortrows([TT']))';
Torg=TTT(2:end,:);
plot(Torg')
xlim([0 c+1]);
xlabel('Monte-Carlo run')
ylabel('\theta')