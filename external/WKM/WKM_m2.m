function [dx, y,varargout] = WKM_m(t, x, u, HbT0,E0,tau,alpha, varargin)
% ODE version of Buxton & Frank's Balloon Model.
% Modified from:
% Modeling the hemodynamic response to brain activation
% Richard B. Buxton,* Kamil Uludag, David J. Dubowitz, and Thomas T. Liu
% NeuroImage 23 (2004) S220 ? S233

% Outputs:
% HbO2
% HbR
%
% States:
% 1: q - normalized HbR
% 2: v - normalized Hb-total
% 3: fin-  flow in
% 4: CMRO2
%
% Inputs:
% 1: delta-Fin
% 2: delta-CMRO2;
%
% Paramaters:
% HbT0 - baseline total-Hb concentration (uM)
% E0   - baseline oxygen extraction fraction
% tau - mean transit through the balloon at rest
% alpha

q=x(1);
v=x(2);
fin=x(3);
CMRO2=x(4);

% fout=(v)^(1/alpha);
% dq_dt = 1/tau*(CMRO2-fout*q/v);
% dv_dt = 1/tau*(Fin-fout);


% 
% q=(CMRO2-tau*dq_dt)*(v)^(1-1/alpha);
% v=(fin-tau*dv_dt)^(alpha);
% fin = tau*dv_dt+(v)^(1/alpha)
% CMRO2=tau*dq_dt+(v)^(1/alpha-1)*q;

% dq_dt = 1/tau*(CMRO2-(v)^(1/alpha-1)*q);
% dv_dt = 1/tau*(Fin-(v)^(1/alpha));

% State equations.
dx(1,1)=1/tau*(CMRO2-(v)^(1/alpha-1)*q); % dq_dt
dx(2,1)=1/tau*(fin-(v)^(1/alpha));  % dq_dt
dx(3,1)=u(1);  %dFin/dt
dx(4,1)=u(2);

% Measurement
y(1,1)=(v-1)*HbT0-(q-1)*HbT0*E0;  %HbO2
y(2,1)=(q-1)*HbT0*E0;  % Hb
  
dF_dx=zeros(4,4);
if(nargout>2)
   % Return dF/dx
  
  % q=(CMRO2-tau*dq_dt)*(v)^(1-1/alpha);
  % q=CMRO2*(v)^(1-1/alpha);
   dF_dx(1,1)= 0;     
   dF_dx(1,2)=x(4)*(1-1/alpha)*(x(2))^(-1/alpha);
   dF_dx(1,3)= 0; 
   dF_dx(1,4)=(x(2))^(1-1/alpha);
   
   %v=(fin-tau*dv_dt)^(alpha);
   %v=(fin)^(alpha);
   dF_dx(2,1)=0; 
   dF_dx(2,2)=0; 
   dF_dx(2,3)=alpha*x(3)^(alpha-1);
   dF_dx(2,4)=0;  
   
   %fin = tau*dv_dt+(v)^(1/alpha)
   %fin = (v)^(1/alpha)
   dF_dx(3,1)=0;
   dF_dx(3,2)=(1/alpha)*x(2)*(1/alpha-1); 
   dF_dx(3,3)=0;
   dF_dx(3,4)=0;
   
   %CMRO2=tau*dq_dt+(v)^(1/alpha-1)*q;
   %CMRO2=(v)^(1/alpha-1)*q;
   dF_dx(4,1)=(x(2))^(1/alpha-1);
   dF_dx(4,2)=(1/alpha-1)*(x(2))^(1/alpha-2)*x(1); 
   dF_dx(4,3)=0;
   dF_dx(4,4)=0;
   
   varargout{1}=dF_dx;
end
  
if(nargin>3)
    %dY_dx  
    %y(1)=(v-1)*HbT0-(q-1)*HbT0*E0;  %HbO2
    %y(2)=(q-1)*HbT0*E0;  % Hb
    
    dY_dx=zeros(2,4);
    dY_dx(1,:)=HbT0*dF_dx(2,:)-dF_dx(1,:)*HbT0*E0;
    dY_dx(2,:)=HbT0*E0*dF_dx(1,:);
    varargout{2}=dY_dx;
end
  
