function nlgr = WKM_idnl(a1,a2,a3,a4,a5,a6)

if(nargin==0)
    a1 = 4; %tau
    a2 = 0.38; %alpha
    a3 = .4; % OEF0
    a4 = 50; %HbT0
end

FileName      = 'WKM_m';   % File describing the WKM model structure.
Order         = [2 1 3];   % Model orders [ny nu nx].
Parameters    = [a1;a2;a3;a4];   % Initial parameters. Np = 6.
InitialStates = [1;1;1];   % Initial initial states.
Ts            = 0;           % Time-continuous system.
nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
    'Name', 'WIndkessel model');

set(nlgr, 'InputName', {'s'}, 'InputUnit', {'%'},...
    'OutputName', {'delta-HbO2', 'delta-HbR'}, ...
    'OutputUnit', {'uM', 'uM'},                         ...
    'TimeUnit', 's');

nlgr = setinit(nlgr, 'Name', {'f','v','q'});
nlgr = setinit(nlgr, 'Unit', {'%','%','%'});
nlgr = setpar(nlgr, 'Name', {'a1','a2','a3','a4'});
nlgr = setpar(nlgr, 'Minimum', {1,.25,.2,30});
nlgr = setpar(nlgr, 'Maximum', {6,.5,.5,70});
nlgr = setpar(nlgr, 'Fixed', {0,0,0,0});

nlgr.SimulationOptions.AbsTol = 1e-6;
nlgr.SimulationOptions.RelTol = 1e-5;

return