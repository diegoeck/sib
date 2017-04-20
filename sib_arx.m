function [theta]=sib_arx(u,y,na,nb,nz)
%  [theta] = sib_arx(u,y,na,nb,nz)


phi=[];

for i=1+nz:nb+nz
    phi=[phi [zeros(i-1,1) ; u(1:end-i+1)]];
    
end

for i=1:na
    phi=[phi [zeros(i,1) ; -y(1:end-i)]];
    
end

theta=(phi'*phi)\(phi'*y);



