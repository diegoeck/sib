function [theta3]=sib_armax(u,y,na,nb,nc,nz)
%  [theta] = sib_oe(u,y,nf,nb,nz)


 [theta] = sib_arx(u,y,na,nb,nz);

 theta=[theta; zeros(nc,1)]
 %theta=[theta; theta(nb+1:nb+na)];
 
 
 %[theta2] = sib_oe_conjugate(u,y,theta,nz,nb); 

 [theta3,J] = sib_armax_hessian(u,y,theta,200,nz,na,nb);




