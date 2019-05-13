function [teta,J] = sib_gradiente(f,teta,N,P)

d=zeros(size(teta));

alfa=1e-5;

[yc,g,J,H] = f(teta,P);

Jv=J;

for i=1:N

    [yc,g,J,H] = f(teta,P);

    d=(4*d+g)/5;

    if J>Jv
        alfa=alfa*0.99;
    else
        alfa=alfa*1.01;
    end
    Jv=J;
    
    passo=0.01/norm(d)*alfa;
    teta=teta-passo*d;

    %disp([i J passo]);

    if mod(i,100)==0
        fprintf('GRADIENTE: J = %d, Alfa = %d, Teta =',J,alfa)
        fprintf('%f ',teta')
        fprintf('\n')

    end
end