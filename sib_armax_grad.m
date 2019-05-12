function [yc,g,J,H] = sib_armax_grad(teta,P)

u=P.u;
y=P.y;
nz=P.nz;
na=P.na;
nb=P.nb;


N=length(u);
M=length(teta);

B=[zeros(nz,1); teta(1:nb)];
A=[1; teta(nb+1:na+nb)];
C=[1; teta(na+nb+1:end)];

yc=filter(B,C,u)+y-filter(A,C,y);
e=yc-y;
J=mean(e.^2);

g=zeros(1,M);
F=[];

for i=1:nb
    f=filter([zeros(nz,1); zeros(i-1,1) ;1 ],C,u);
    F=[F f] ;
end

for i=1:na
    f=-filter([0; zeros(i-1,1) ;1 ],C,y);
    F=[F f] ;
end


for i=1:(M-na-nb)
    f=-filter(B,C,u)+filter(A,C,y);
    f=filter([0; zeros(i-1,1) ;1 ],C,f);
    F=[F f] ;
end

g=(F'*e)/N;

H=zeros(M,M);

for j=1:N
    H=H+F(j,:)'*F(j,:);
end

H=H/N;

end
