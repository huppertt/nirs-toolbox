function nlgr = WKM_idnl(b0,b1,a0,a1,a2,a3,a4,a5,a6)

if(nargin==0)
    b0=1;
    b1=2;
    a0 = 1;  % Gain on flow inducing signal
    a1 = 3; % 1/taus
    a2 = 2; % 1/tauf
    a3 = 4; %tau
    a4 = 0.38; %alpha
    a5 = .4; % OEF0
    a6 = 100; %HbT0
end

FileName      = 'WKM_vs2_m';   % File describing the WKM model structure.
Order         = [2 2 7];   % Model orders [ny nu nx].
Parameters    = [b0;b1; a0;a1;a2;a3;a4;a5;a6];   % Initial parameters. Np = 7.
InitialStates = [0;0;1;1;1;1;1];   % Initial initial states.
Ts            = 0;           % Time-continuous system.
nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
    'Name', 'Windkessel model');

set(nlgr, 'InputName', {'flow-inducing','CMRO2-inducing'}, 'InputUnit', {'%','%'},...
    'OutputName', {'delta-HbO2', 'delta-HbR'}, ...
    'OutputUnit', {'uM', 'uM'},                         ...
    'TimeUnit', 's');

nlgr = setinit(nlgr, 'Name', {'flow-inducing','CMRO2-inducing','CBF','CMRO2','CBV','OEF','q'});
nlgr = setinit(nlgr, 'Unit', {'%','%','%','%','%','%','%'});
nlgr = setpar(nlgr, 'Name', {'b0','b1','a0','a1','a2','a3','a4','a5','a6'});
nlgr = setpar(nlgr, 'Minimum', {-Inf,1,-Inf,1,1,1,.25,.2,1});
nlgr = setpar(nlgr, 'Maximum', {Inf,6,Inf,6,6,6,.5,.7,200});
nlgr = setpar(nlgr, 'Fixed', {1,1,1,1,1,1,1,1,1});

nlgr.SimulationOptions.AbsTol = 1e-6;
nlgr.SimulationOptions.RelTol = 1e-5;

return