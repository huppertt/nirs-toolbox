classdef  WKM_static < nirs.vascular.models.AbstractVascularModel
   
    
    methods
        function obj = WKM_static(varargin)
            obj.model=@nirs.vascular.models.WKM.static;
            obj.description = 'Steady-state vascular model'
            obj.fitfunction = 'static';
            obj.states={'CMRO2','CBF','OEF','CBV','q'};
            
            obj.parameters.tau=3;
            obj.parameters.alpha=.38;
            obj.parameters.HbT0=100;
            obj.parameters.E0=.40;
            
            obj.outputNames={'CMRO2','CBF'};
             
        end
        
        
    end
end