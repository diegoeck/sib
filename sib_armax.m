function [theta3]=sib_oe(u,y,na,nb,nz)
%  [theta] = sib_oe(u,y,nf,nb,nz)


 [theta] = sib_arx(u,y,na,nb,nz);

 %[theta2] = sib_oe_conjugate(u,y,theta,nz,nb); 

 [theta3,J] = sib_armax_hessian(u,y,theta,2000,nz,nb);




