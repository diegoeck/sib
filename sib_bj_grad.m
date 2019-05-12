function [yc,g,J,H] = sib_armax_grad(teta,P)

u=P.u;
y=P.y;
nz=P.nz;
na=P.na;
nb=P.nb;
nc=P.nc;

N=length(u);
M=length(teta);

B=[zeros(nz,1); teta(1:nb)];
A=[1; teta(nb+1:na+nb)];
C=[1; teta(na+nb+1:na+nb+nc)];

D=[1; teta(na+nb+nc+1:end)];

yc=filter(B,A,u)-y;
yc=filter(D,C,yc)+y;

e=yc-y;
J=mean(e.^2);

g=zeros(1,M);
F=[];
    
uf=filter(D,C,u);
for i=1:nb
    f=filter([zeros(nz,1); zeros(i-1,1) ;1 ],A,uf);
    F=[F f] ;
end

for i=1:na
    yf=filter(B,A,uf);
    f=-filter([0; zeros(i-1,1) ;1 ],A,yf);
    F=[F f] ;
end

yf=filter(B,A,u)-y;
for i=1:nc
    f=filter(D,C,yf);
    f=-filter([0; zeros(i-1,1) ;1 ],C,f);
    F=[F f] ;
end

for i=1:(M-na-nb-nc)
    f=filter([0; zeros(i-1,1) ;1 ],C,yf);
    F=[F f] ;
end


g=(F'*e)/N;

H=zeros(M,M);

for j=1:N
    H=H+F(j,:)'*F(j,:);
end

H=H/N;

end
