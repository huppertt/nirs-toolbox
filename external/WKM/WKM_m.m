function [dx, y,varargout] = WKM_m(t, x, u, a1,a2,a3,a4, varargin)
% State space model of the hemodynamic balloon model.  Modified from Riera
% et al NeuroImage 21 (2004) 547-567.

% States
% x1 - Flow-in
% x2 - v
% x3 - q

%a1 = tau = 4s
%a2 = alpha = 0.38
%a3 = OEF0 = .40 
%a4 = HbT0 = 50uM

dx(1,1) = u(1);  % Flow
dx(2,1) = 1/a1*(x(1)-x(2)^(1/a2));  % CBV
dx(3,1) = 1/a1*(x(1)/a3*(1-(1-a3)^(1/x(2)))-x(3)*x(2)^((1-a2)/a2));


% Measurement
y(1,1)=(x(2)-1)*a4-(x(3)-1)*a3*a4;  %HbO2
y(2,1)=(x(3)-1)*a3*a4;  % Hb
  
Jx=zeros(4,4);
if(nargout>2)
    % Return dF/dx
%     
%     Jx(1,1)=-a1;
%     Jx(1,2)=-a2;
%     Jx(2,1)=1;
%     Jx(3,2)=1/a3;
%     Jx(3,3)=-x(3)^((1-a4)/a4)/(a3*a4);
%     Jx(4,2)=1/(a3*a5)*(1-(1-a5)^(1/x(2))+log(1-a5)*(1-a5)^(1/x(2))/x(2));
%     Jx(4,3)=-(1-a4)*x(4)*x(3)^((1-2*a4)/a4)/(a3*a4);
%     Jx(4,4)=-x(3)^((1-a4)/a4)/a3;
%     
    varargout{1}=Jx;
end
dY_dx=zeros(2,4);
 
if(nargin>3)
    %dY_dx  
    %y(1)=(v-1)*HbT0-(q-1)*HbT0*E0;  %HbO2
    %y(2)=(q-1)*HbT0*E0;  % Hb
    
    dY_dx(1,:)=a4*Jx(2,:)-Jx(3,:)*a3*a4;
    dY_dx(2,:)=Jx(2,:)*a3*a4;
    varargout{2}=dY_dx;
end
  
