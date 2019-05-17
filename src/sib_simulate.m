function [ys] = sib_simulate(u, m)
% [ys] = sib_simulate(u, m)
%
% Simulate model
%    
% ys(t) = G(z)*u(t)

ys = filter(m.B, m.A, u);
ys = filter(1, m.F, ys);


