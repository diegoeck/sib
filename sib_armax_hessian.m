function [teta,J] = sib_armax_hessian(u,y,teta,NN,nz,na,nb)

d=zeros(size(teta));

for i=1:10

    [yc,g,J,H] = sib_armax_grad(u,y,teta,nz,na,nb);

    %d=(4*d+g)/5;
    
    passo=0.01/norm(g);
    teta=teta-passo*g;
    
    %disp([i J passo]);  
    
    if mod(i,100)==0
        sprintf('J %f %f',J,teta')
    end
end

for i=1:NN

    [yc,g,J,H] = sib_armax_grad(u,y,teta,nz,na,nb);

    passo=i/NN;
    teta=teta-passo*(H\g);

    if mod(i,100)==0
        sprintf('J %f %f',J,teta')
    end
    
end
