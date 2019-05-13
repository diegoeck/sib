function [theta, m, estavel] = sib_oe_filtered(u, y, na, nb, nz);
%  [theta, m] = sib_oe_filtered(u, y, nf, nb, nz)
%
%  Prediciton error method with OE structure
%
%         B(z)          
%  y(t) = ---- u(t) + e(t) 
%         F(z)        
%

FF=[.1 .2 .3 .4 .5 .6 .7 .8 .9];

[thetai] = sib_arx(u,y,na,nb,nz);

theta=thetai;

for i=1:length(FF)

    %Filter the data
    [fb,fa] = butter(1,FF(i));
    yf=filter(fb,fa,y);
    uf=filter(fb,fa,u);  

    %Estimate with filtered data

    [theta2] = sib_oe_c(uf,yf,theta,nb,nz) ;
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
[thetaf] = sib_oe_c(u,y,theta,nb,nz) ;

m.A=1;
m.B=[zeros(nz,1); theta2(1:nb)];
m.C=1;
m.D=1;
m.F=[1; theta2(nb+1:end)];


