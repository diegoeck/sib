function [theta]=sib_sm(u,y,na,nb,nz)

[theta] = sib_arx(u,y,na,nb,nz);


for i=1:10
%F=tf([1 zeros(1,na)],[1; theta(nb+1:nb+na)]',1);
uf=filter(1,[1; theta(nb+1:nb+na)]',u);
yf=filter(1,[1; theta(nb+1:nb+na)]',y);

    
[theta] = sib_arx(uf,yf,na,nb,nz);
    
    
end


