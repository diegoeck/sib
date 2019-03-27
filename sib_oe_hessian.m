function [teta,J] = sib_oe_hessian(u,y,teta,NN,nz,nb)

for i=1:NN

    [yc,g,J,H] = sib_oe_grad(u,y,teta,nz,nb);

    passo=i/NN;
    
    q=H\g;
    
    teta=teta-passo*q;
    
    if mod(i,10)==0
        k=floor(i/100);
        fprintf('\r %1.10f %1.10f ',J,passo)
        for j=1:k
            fprintf('#')
        end
        for j=1:(20-k)
            fprintf('-')
        end
    end
    
end

