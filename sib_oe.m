function [theta2]=sib_oe(u,y,na,nb,nz)
%  [theta] = sib_oe(u,y,nf,nb,nz)

[theta] = sib_arx(u,y,na,nb,nz);
[theta2] = sib_oe_c(u,y,theta,nb,nz)





