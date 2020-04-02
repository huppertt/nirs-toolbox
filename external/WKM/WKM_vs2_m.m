function [dx, y,varargout] = WKM_m(t, x, u, b0,b1,a0,a1,a2,a3,a4,a5,a6, varargin)
% State space model of the hemodynamic balloon model.  Modified from 
% Neuroimage. 2004 Feb;21(2):547-67.
% A state-space model of the hemodynamic approach: nonlinear filtering of BOLD signals.
% Riera JJ1, Watanabe J, Kazuki I, Naoki M, Aubert E, Ozaki T, Kawashima R.

% States
% x1 - flow indcing signal
% x2 - metabolism indcing signal
% x3 - f - CBF(t)/CBF0
% x4 - m - CMRO2(t)/CMRO20
% x5 - v - CBV(t)/CBV0
% x6 - e - OEF(t)/OEF0
% x7 - q - HbR(t)/HbR0
%b0 = gain on cmro2 inducing signal
%b1 = taus = 4s
%a0 = gain on flow inducing signal
%a1 = taus = 4s
%a2 = tauf = 4s - autoregulation feedback
%a3 = tau = 4s - 
%a4 = alpha = 0.38
%a5 = OEF0 = .40 
%a6 = HbT0 = 50uM

% MATHEMATICA NOTATION:
% F1 = a0*u1 - a1*x1 - a2*(x3 - 1)
% F2 = b0*u2 - b1*x2 - a2*(x4 - 1)
% F3 = x1
% F4 = x2
% F5 = 1/a3*(x3 - x5^(1/a4))
% F6 = (x3*x2 - x4*x1)/x3^2
% F7 = 1/a3*x3*(x6 - x7/x5) + x7/x5*(1/a3*(x3 - x5^(1/a4)))
% F = {F1, F2, F3, F4, F5, F6, F7}


% Equations 2-5
dx(1,1) = a0*u(1) - a1*x(1)-a2*(x(3)-1);  % flow inducing signal (w/ autoregulation feedback )
dx(2,1) = b0*u(2) - b1*x(2)-a2*(x(4)-1);  % cmro2 inducing signal

dx(3,1) = x(1);  % Flow
dx(4,1) = x(2);  % CMRO2

dx(5,1) = 1/a3*(x(3)-x(5)^(1/a4));  % CBV

% e = m/f
dx(6,1) = (x(3)*x(2)-x(4)*x(1))/x(3)^2;    
       
% dq/dt = 1/tau*f*(e-q/v)+q/v*dv/dt
dx(7,1) = 1/a3*x(3)*(x(6)-x(7)/x(5))+x(7)/x(5)*dx(5);

% Measurement
y(1,1)=(x(5)-1)*a6-(x(7)-1)*a5*a6;  %HbO2
y(2,1)=(x(7)-1)*a5*a6;  % Hb
 

Jx=zeros(7,7);
if(nargout>2)
    % Return dF/dx
    
% MATHEMATICA NOTATION:
% D[F, {{x1, x2, x3, x4, x5, x6, x7}}]
    
    Jx(1,1)=-a1;
    Jx(1,3)=-a2;
    Jx(2,2)=-b1;
    Jx(2,4)=-a2;
    Jx(3,1)=1;
    Jx(4,2)=1;
    Jx(5,3)=1/a3;
    Jx(5,5)=-(x(5)^(-1+1/a4)/(a3*a4));
    Jx(6,1)=-(x(4)/x(3)^2);
    Jx(6,2)=1/x(3);
    Jx(6,3)=x(2)/x(3)^2-(2*(x(2)*x(3)-x(1)*x(4)))/x(3)^3;
    Jx(6,4)=-(x(1)/x(3)^2);
    Jx(7,3)=x(7)/(a3*x(5))+(x(6)-x(7)/x(5))/a3;
    Jx(7,5)=(x(3)*x(7))/(a3*x(5)^2)-(x(5)^(-2+1/a4)*x(7))/(a3*a4)-((x(3)-x(5)^(1/a4))*x(7))/(a3*x(5)^2);
    Jx(7,6)=x(3)/a3;
    Jx(7,7)=-(x(3)/(a3*x(5)))+(x(3)-x(5)^(1/a4))/(a3*x(5));
    Jx=Jx+eye(7);
    varargout{1}=Jx;
end
dY_dx=zeros(2,7);
 
if(nargout>3)
    % MATHEMATICA NOTATION:
    % Y1 = (x5 - 1)*a6 - (x7 - 1)*a5*a6
    % Y2 = (x7 - 1)*a5*a6
    % Y = {{Y1, Y2}}
    % D[Y, {{x1, x2, x3, x4, x5, x6, x7}}]
    
    %dY_dx  
    %y(1)=(v-1)*HbT0-(q-1)*HbT0*E0;  %HbO2
    %y(2)=(q-1)*HbT0*E0;  % Hb
    
    dY_dx=[0 0 0 0 a6 0 -a5*a6;...
           0 0 0 0 0  0  a5*a6]*Jx;
    varargout{2}=dY_dx;
end
  
