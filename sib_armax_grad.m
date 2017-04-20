function [g,H,J] = sib_armax_grad(u,y,teta,nz,nb,nc)

N=length(u);
M=length(teta);

B=[zeros(nz,1); teta(1:nb)];
A=[1; teta(nb+1:M-nc)];
C=[1; teta(M-nc+1:end)];

%yp=filter(B,C,u)+y-filter(A,C,y) 

e=+filter(B,C,u)-filter(A,C,y) ;
J=mean(e.^2);

g=zeros(1,M);
F=[];

fu=filter(1,C,u);
for i=1:nb
    F=[F [zeros(nz,1); zeros(i-1,1);fu(1:end-nz-i+1)]]; 
end

fu=filter(1,C,y);
for i=1:M-nb-nc
    F=[F [0; zeros(i-1,1);-fu(1:end-i)]]; 
end


fu=filter(B,C,u);
fy=filter(A,C,y);
fy=filter(1,C,fu-fy);
for i=1:nc
    F=[F [0; zeros(i-1,1); -fy(1:end-i)]];
end

g=(F'*e)/N;

H=zeros(M,M);

for j=1:N
    H=H+F(j,:)'*F(j,:);
end

H=H/N;

end