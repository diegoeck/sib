function [teta,J] = sib_armax_hessian(u,y,teta,NN,nz,na,nb)

d=zeros(size(teta));

alfa=1e-5;

    [yc,g,J,H] = sib_armax_grad(u,y,teta,nz,na,nb);

    Jv=J;
for i=1:10000

    [yc,g,J,H] = sib_armax_grad(u,y,teta,nz,na,nb);

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
        disp([J teta']);
        disp(alfa);
        %fprintf('J %f teta %f \n',J,teta')
    end
end

for i=1:NN

    [yc,g,J,H] = sib_armax_grad(u,y,teta,nz,na,nb);

    passo=i/NN;
    teta=teta-passo*(H\g);

    if mod(i,100)==0
        disp([J teta']);

        %fprintf('J %f teta %f \n',J,teta')
    end

end
