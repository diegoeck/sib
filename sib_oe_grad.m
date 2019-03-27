function [yc,g,J,H] = sib_oe_grad(u,y,teta,nz,nb)

N=length(u);
M=length(teta);

B=[zeros(nz,1); teta(1:nb)];
A=[1; teta(nb+1:end)];

yc=filter(B,A,u); 
e=yc-y;
J=mean(e.^2);

g=zeros(1,M);
F=[];


f=filter(1,A,u);
for i=1:nb
    F=[F [zeros(nz,1); zeros(i-1,1);f(1:end-nz-i+1)]]; 
end




%f=filter([0; zeros(i-1,1) ;1 ],A,yc);
f=filter([1 ],A,yc);

for i=1:(length(teta)-nb)
    F=[F [0; zeros(i-1,1); -f(1:end-i)] ];
end

g=(F'*e)/N;

H=zeros(M,M);

for j=1:N
    H=H+F(j,:)'*F(j,:);
end

H=H/N;

end