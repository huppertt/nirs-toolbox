classdef  WKM_static < nirs.vascular.models.AbstractVascularModel
    methods
        function obj = WKM_static(varargin)
            obj.model=@nirs.vascular.models.WKM.static;
                         
            obj.parameters(1).name='tau';
            obj.parameters(1).value= 3;
            obj.parameters(1).fixed= false;
            obj.parameters(1).range= [1 5];
            
            obj.parameters(2).name='alpha';
            obj.parameters(2).value= .38;
            obj.parameters(2).fixed= false;
            obj.parameters(2).range= [.2 .5];
            
            obj.parameters(3).name='HbT0';
            obj.parameters(3).value= 100;
            obj.parameters(3).fixed= false;
            obj.parameters(3).range= [30 200];
            
            obj.parameters(4).name='E0';
            obj.parameters(4).value= .40;
            obj.parameters(4).fixed= false;
            obj.parameters(4).range= [.2 .65];
            
            obj.name='Static Model';
            obj.description = 'Steady-state vascular model'
            obj.outputNames= {'CMRO2','CBF'};  % Default output names
            obj.states={'CBF','OEF','CMRO2','CBV','q'};
            
            obj.model=model(obj);
            
        end
        
        
       function nlgr= model(obj)
            
            paramorder={'tau','alpha','HbT0','E0'};
            
            for i=1:length(paramorder)
                id=find(ismember({obj.parameters.name},paramorder{i}));
                pValues(i)=obj.parameters(id).value;
                pValuesMin{i}=obj.parameters(id).range(1);
                pValuesMax{i}=obj.parameters(id).range(2);
                pValuesFixed{i}=obj.parameters(id).fixed;
            end
            
            FileName      = @nirs.vascular.models.WKM.WKM_static;   % File describing the WKM model structure.
            Order         = [2 2 4];   % Model orders [ny nu nx].
            Parameters    = pValues;   % Initial parameters. Np = 7.
            InitialStates = [0; 0;1;1;1;1;1];   % Initial initial states.
            Ts            = 0;           % Time-continuous system.
            nlgr = idnlgrey(FileName, Order, Parameters, InitialStates, Ts, ...
                'Name', 'Windkessel model');
            
            set(nlgr, 'InputName', {'flow-inducing','CMRO2-inducing'}, 'InputUnit', {'%','%'},...
                'OutputName', {'HbO2', 'HbR'}, ...
                'OutputUnit', {'uM', 'uM'},                         ...
                'TimeUnit', 's');
            
            nlgr = setinit(nlgr,'Name', obj.states);
            nlgr = setpar(nlgr, 'Name', paramorder);
            nlgr = setpar(nlgr, 'Minimum', pValuesMin);
            nlgr = setpar(nlgr, 'Maximum', pValuesMax);
            nlgr = setpar(nlgr, 'Fixed', pValuesFixed);
            
            nlgr.SimulationOptions.AbsTol = 1e-6;
            nlgr.SimulationOptions.RelTol = 1e-5;
            
        end
        
        
    end
end