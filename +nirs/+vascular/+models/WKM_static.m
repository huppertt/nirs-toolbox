classdef  WKM_static < nirs.vascular.models.AbstractVascularModel
    methods
        function obj = WKM_static(varargin)
            obj.model=@nirs.vascular.models.WKM.static;
             
            obj.parameters.tau=3;
            obj.parameters.alpha=.38
            obj.parameters.HbT0=100;
            obj.parameters.E0=.40;
            
            obj.name='Static Model';
            obj.description = 'Steady-state vascular model'
            obj.outputNames= {'CMRO2','CBF'};  % Default output names
            obj.states={'CBF','OEF','CMRO2','CBV','q'};
            
            
        end
        
    end
end