function [f] = ksi_calc_sp1(mesh)

%Calculates the boundary condition values for SP3 case based on equations 28a & 28b
%of Klose 2006 paper. 

BC=bound_coeffs(mesh,1);

A=BC(:,1);
B=BC(:,2);

f = (1+B(:,1)).\((1/2)+A(:,1));
