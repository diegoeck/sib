function teta=oe(u,y,teta,NN,nz,nb)

for i=1:NN
    [yc,g,J,H] = sib_oe_grad(u,y,teta,nz,nb);
    passo=0.001/norm(g);
    teta=teta+passo*g;
end

for i=1:100
    [yc,g,J,H] = sib_oe_grad(u,y,teta,nz,nb);
    teta=teta+H\g/(101-i);
end


for i=1:100
    [yc,g,J,H] = sib_oe_grad(u,y,teta,nz,nb);
    teta=teta+H\g;
end



end