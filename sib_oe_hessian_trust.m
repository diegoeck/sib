function [teta,J] = sib_oe_hessian_trust(u,y,teta,NN,nz,nb)

passo=1e-6;

[yc,g,J,H] = sib_oe_grad(u,y,teta,nz,nb);

for i=1:NN

    [yc,g,Jn,H] = sib_oe_grad(u,y,teta,nz,nb);

    if Jn>J
        passo=passo*0.99;
    else
        passo=passo*1.01;
    end
    
    
    J=Jn;
        
    dir=H\g;
    dir=dir/norm(dir);
    dir=dir*passo;
    
    teta=teta+dir;
    
    disp([i J passo])  
    
    if passo>1
        break
    end
    
end

for i=1:50

    [yc,g,Jn,H] = sib_oe_grad(u,y,teta,nz,nb);

        
    dir=H\g;
    
    teta=teta+dir;
    
    %disp([i J passo])  
    
    
end
