function [teta,J] = sib_armax_hessian(u,y,teta,NN,nz,nb,nc)

for i=1:NN

    [g,H,J] = sib_armax_grad(u,y,teta,nz,nb,nc);

    passo=i/NN;
        
    teta=teta-passo*(H\g)
    %teta=teta-(H\g)/1000
    
    
    %disp([i J passo]) 
    %disp(J)
    
end
