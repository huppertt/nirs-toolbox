function [f1,f2,g1,g2] = ksi_calc_sp3(mesh)

%Calculates the boundary condition values for SP3 case based on equations 28a & 28b
%of Klose 2006 paper.  

[BC]=bound_coeffs(mesh,3);

A=BC(:,1:2);
B=BC(:,3:4);
C=BC(:,5:6);
D=BC(:,7:8);

X1=((7/24)+A(:,2))-3.*D(:,2).*((1+B(:,1)).\((1/8)+C(:,1)));
X2=(1+B(:,2))-21.*D(:,2).*((1+B(:,1)).\D(:,1));
X3=((1/8)+C(:,2))-3.*D(:,2).*((1+B(:,1)).\((1/2)+A(:,1)));

f1=(1+B(:,1)).\(((1/2)+A(:,1))-7.*(D(:,1).*(X2.\X3)));
f2=X2.\X1;
g1=-(1+B(:,1)).\(((1/8)+C(:,1))-7.*D(:,1).*(X2.\X1));
g2=-X2.\X3;


