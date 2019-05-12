function [teta,J] = sib_newton(f,teta,N,P)

for i=1:N

    [yc,g,J,H] = f(teta,P);

    passo=i/N;
    teta=teta-passo*(H\g);

    if mod(i,100)==0
        fprintf('NEWTON: J = %d, Teta =',J)
        fprintf('%f ',teta')
        fprintf('\n')
    end

end

