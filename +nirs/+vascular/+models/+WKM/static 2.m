function [states,yhat] = static(t,HbO,HbR)
% Steady-state version of WKM

tau=param.tau;
alpha=param.alpha;
HbT0=param.HbT0;
E0=param.E0;


q = 1+HbR/(E0*HbT0);
v = 1+(HbO+HbR)/HbT0;

dq = [0; diff(q)./diff(t)];
dv = [0; diff(v)./diff(t)];

f = v.^(1/alpha);
OEF = q./v+tau./f.*(dq-q./v.*dv);
m = f.*OEF;

states(1).name='CBF';
states(1).data=f;
states(2).name='OEF';
states(2).data=OEF;
states(3).name='CMRO2';
states(3).data=m;
states(4).name='CBV';
states(4).data=v;
states(5).name='q';
states(5).data=q;

yhat=[HbO2 HbR];

return