function [theta,estavel]=sib_oe_filtered(u,y,na,nb,nz);
%

FF=[.1];

[thetai] = sib_arx(u,y,na,nb,nz);

theta=thetai;


for i=1:length(FF)

    %Filter the data
    [fb,fa] = butter(1,FF(i));
    yf=filter(fb,fa,y);
    uf=filter(fb,fa,u);  

    %Estimate with filtered data

    [theta2] = sib_oe_conjugate(uf,yf,theta,nz,nb) ;
%    [thetaf,J] = sib_oe_hessian(uf,yf,theta2,2000,nz,nb);
thetaf=theta2;

    %Test if the model is stable
    if isnan(thetaf(1))
        theta=thetai;
        estavel(i)=1;
        
    elseif max(abs(roots([1;thetaf(nb+1:na+nb)])))>1
        theta=thetai;
        estavel(i)=1;
    else
        theta=thetaf;
        estavel(i)=0;
    end    
    
end


%Estimate with real data
[thetaf] = sib_oe_conjugate(u,y,theta,nz,nb) ;
[theta,J] = sib_oe_hessian(u,y,thetaf,50,nz,nb);

%theta=thetaf;

