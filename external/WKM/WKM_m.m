function [dx, y,varargout] = WKM_m(t, x, u, a0,a1,a2,a3,a4,a5,a6, varargin)
% State space model of the hemodynamic balloon model.  Modified from 
% Neuroimage. 2004 Feb;21(2):547-67.
% A state-space model of the hemodynamic approach: nonlinear filtering of BOLD signals.
% Riera JJ1, Watanabe J, Kazuki I, Naoki M, Aubert E, Ozaki T, Kawashima R.


% States
% x1 - Flow-in
% x2 - v
% x3 - q

%a1 = taus = 4s
%a2 = tauf = 4s
%a3 = tau = 4s
%a4 = alpha = 0.38
%a5 = OEF0 = .40 
%a6 = HbT0 = 50uM

% Equations 2-5  (note a1 is my 1/a1)
dx(1,1) = a0*u(1) - a1*x(1)-a2*(x(2)-1);  % flow inducing signal (w/ feedback shutoff)
dx(2,1) = x(1);  % Flow
dx(3,1) = 1/a3*(x(2)-x(3)^(1/a4));  % CBV
dx(4,1) = 1/a3*(x(2)/a5*(1-(1-a5)^(1/x(2)))-x(4)*x(3)^((1-a5)/a5));


% Measurement
y(1,1)=(x(3)-1)*a6-(x(4)-1)*a5*a6;  %HbO2
y(2,1)=(x(4)-1)*a5*a6;  % Hb
  
Jx=zeros(4,4);
if(nargout>2)
    % Return dF/dx
    
    Jx(1,1)=-a1;
    Jx(1,2)=-a2;
    Jx(2,1)=1;
    Jx(3,2)=1/a3;
    Jx(3,3)=-x(3)^((1-a4)/a4)/(a3*a4);
    Jx(4,2)=1/(a3*a5)*(1-(1-a5)^(1/x(2))+log(1-a5)*(1-a5)^(1/x(2))/x(2));
    Jx(4,3)=-(1-a4)*x(4)*x(3)^((1-2*a4)/a4)/(a3*a4);
    Jx(4,4)=-x(3)^((1-a4)/a4)/a3;
%     
    Jx=Jx+eye(4);
    varargout{1}=Jx;
end
dY_dx=zeros(2,4);
 
if(nargout>3)
    %dY_dx  
    %y(1)=(v-1)*HbT0-(q-1)*HbT0*E0;  %HbO2
    %y(2)=(q-1)*HbT0*E0;  % Hb
    
    dY_dx(1,:)=Jx(3,:)*a6-Jx(4,:)*a5*a6;
    dY_dx(2,:)=Jx(4,:)*a5*a6;
    varargout{2}=dY_dx;
end
  
