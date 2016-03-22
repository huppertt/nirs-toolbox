classdef  WKM_difflimit < nirs.vascular.models.AbstractVascularModel
    
    properties(Access=protected)
        nlgr;
    end
    methods
        function obj = WKM_difflimit(varargin)
            
            obj.description = 'O2 diffusion limited model'
            obj.fitfunction = 'kalman';
            obj.states={'Flow Inducing','CMRO2','CBF','OEF','CBV','q'};
            
            obj.parameters.tau_autoreg = 3;
            obj.parameters.tau_flowind = 2;
            obj.parameters.gain_flowind = 1;
            obj.parameters.tau=3;
            obj.parameters.alpha=.38;
            obj.parameters.HbT0=100;
            obj.parameters.E0=.40;
            
            obj.outputNames={'CMRO2','CBF'};
            obj.model=obj.nlgr;
        end
    
    function nlgr= get.nlgr(obj)
            nlgr =  nirs.vascular.models.WKM.WKM_idnl(obj.parameters);
        end
    end
end