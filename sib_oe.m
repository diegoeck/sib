function [theta3]=sib_oe(u,y,na,nb,nz)
%  [theta] = sib_oe(u,y,nf,nb,nz)


 [theta] = sib_arx(u,y,na,nb,nz);

 [theta2] = sib_oe_gradient(u,y,theta,nz,nb); 

 [theta3,J] = sib_oe_hessian(u,y,theta2,2000,nz,nb);




