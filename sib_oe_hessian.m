function [teta,J] = sib_oe_hessian(u,y,teta,NN,nz,nb)

for i=1:NN

    [yc,g,J,H] = sib_oe_grad(u,y,teta,nz,nb);

    passo=i/NN;
        
    teta=teta+passo*inv(H)*g;
    
    disp([i J passo])    
    
end
