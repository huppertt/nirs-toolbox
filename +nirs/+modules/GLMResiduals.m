classdef GLMResiduals < nirs.modules.AbstractModule
%% GLMResiduals - this extracts the residuals from GLM
%
% Options:
%     GLMjob     - a GLM job (e.g., GLMjob = nirs.modules.GLM)
% Example:
%     j = nirs.modules.GLM();
%     j.type = 'AR-IRLS';
%     b = Dictionary();
%     b('default') = nirs.design.basis.Canonical(); % default basis
%     b('A')       = nirs.design.basis.Gamma();     % a different basis for condition 'A'
%     
%     j.basis = b;
%     
%     j.trend_func = @(t) nirs.design.trend.legendre(t, 3); % 3rd order polynomial for trend
%
%     stats = j.run( hb );
%     
%     job_resid = nirs.modules.GLMResiduals();
%     job_resid.GLMjob = j;
%     
%     hb_residuals = job_resid.run( hb );
%     

    properties
        GLMjob;
    end
    methods
        function obj = GLMResiduals( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM model residuals';
            obj.GLMjob = nirs.modules.GLM();
        end
       
        function data = runThis( obj, data )
           
            S = obj.GLMjob.run( data );
            
            for i = 1:numel(data)
                
                t       = data(i).time;
                stims   = data(i).stimulus;
                nchan   = height(unique(data(i).probe.link,'rows'));
                ncond   = length(S(i).beta)/nchan;
            
                X = nirs.design.createDesignMatrix( stims, t, obj.GLMjob.basis );

                data(i).data = data(i).data - X * reshape(S(i).beta,[nchan ncond])';
                
            end
            
        end
        
    end
    
end

