function nlgr = WKM_idnl(param)

a0=param.gain_flowind;
a1=param.tau_autoreg;
a2 = param.tau_flowind; % 1/tauf
a3 = param.tau; %tau
a4 = param.alpha; %alpha
a5 = param.E0; % OEF0
a6 = param.HbT0; %HbT0

FileName      = 'WKM_m';   % File describing the WKM model structure.
Order         = [2 1 4];   % Model orders [ny nu nx].
Parameters    = [a0; a1;a2;a3;a4;a5;a6];   % Initial parameters. Np = 7.
InitialStates = [0; 1;1;1];   % Initial initial states.
Ts            = 0;           % Time-continuous system.
nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
    'Name', 'WIndkessel model');

set(nlgr, 'InputName', {'flow-inducing'}, 'InputUnit', {'%'},...
    'OutputName', {'HbO2', 'HbR'}, ...
    'OutputUnit', {'uM', 'uM'},                         ...
    'TimeUnit', 's');

nlgr = setinit(nlgr, 'Name', {'flow-inducing','CBF','CBV','q'});
nlgr = setinit(nlgr, 'Unit', {'%','%','%','%'});
nlgr = setpar(nlgr, 'Name', {'a0','a1','a2','a3','a4','a5','a6'});
nlgr = setpar(nlgr, 'Minimum', {-Inf,1,1,1,.25,.2,1});
nlgr = setpar(nlgr, 'Maximum', {Inf,6,6,6,.5,.7,200});
nlgr = setpar(nlgr, 'Fixed', {1,1,1,1,1,1,1});

nlgr.SimulationOptions.AbsTol = 1e-6;
nlgr.SimulationOptions.RelTol = 1e-5;

return