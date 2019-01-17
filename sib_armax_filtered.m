function [theta3]=sib_armax_filtered(u,y,na,nb,nc,nz);
%

FF=[.1 .2 .3 .4 .5 .6 .7 .8 .9];

[theta] = sib_arx(u,y,na,nb,nz);
theta=[theta; zeros(nc,1)]




for i=1:length(FF)

    %Filter the data
    [fb,fa] = butter(1,FF(i));
    yf=filter(fb,fa,y);
    uf=filter(fb,fa,u);  

    %Estimate with filtered data

    [theta3,J] = sib_armax_hessian(uf,yf,theta,200,nz,na,nb);

    theta=theta3;

    
end


%Estimate with real data
[theta3,J] = sib_armax_hessian(uf,yf,theta,200,nz,na,nb);

%theta=thetaf;

