function [yc,g,J,H] = sib_oe_grad(u,y,teta,nz,nb)

N=length(u);
M=length(teta);

B=[zeros(nz,1); teta(1:nb)];
A=[1; teta(nb+1:end)];

yc=filter(B,A,u); 
e=y-yc;
J=mean(e.^2);

g=zeros(1,M);
F=[];

for i=1:nb
    f=filter([zeros(nz,1); zeros(i-1,1) ;1 ],A,u); 
    F=[F f] ;  
end

for i=1:(length(teta)-nb)
    f=filter([0; zeros(i-1,1) ;1 ],A,yc); 
    F=[F -f] ;  
end

g=(F'*e)/N;

H=zeros(M,M);

for j=1:N
    H=H+F(j,:)'*F(j,:);
end

H=H/N;

end