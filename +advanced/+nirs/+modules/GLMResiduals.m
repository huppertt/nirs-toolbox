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
        remove_events;
    end
    methods
        function obj = GLMResiduals( prevJob )
            if nargin > 0, obj.prevJob = prevJob; end
            
            obj.name = 'GLM model residuals';
            obj.GLMjob = nirs.modules.GLM();
            obj.remove_events=false;
        end
       
        function data = runThis( obj, data )
           
           
            if(obj.remove_events)
                dorig=data;
                for i = 1:numel(data)
                    job=nirs.modules.KeepStims;
                    keys=nirs.getStimNames(data(i));
                    for k=1:length(keys)
                        if(isa(data(i).stimulus(keys{k}),'nirs.design.StimulusVector'));
                            job.listOfStims{end+1,1}=keys{k};
                        end
                    end
                    if(length(job.listOfStims)==0)
                        warning('cannot remove all events');
                    end
                    data(i)=job.run(data(i));
                end
            end
            
            S = obj.GLMjob.run( data );
            
            for i = 1:numel(data)
                
                t       = data(i).time;
                stims   = data(i).stimulus;
                nchan   = height(unique(data(i).probe.link,'rows'));
                ncond   = length(S(i).beta)/nchan;
            
                X = nirs.design.createDesignMatrix( stims, t, obj.GLMjob.basis );

                data(i).data = data(i).data - X * reshape(S(i).beta,[nchan ncond])';
                
                
                if(obj.remove_events)
                    data(i).stimulus=dorig(i).stimulus;
                    
                    job=nirs.modules.KeepStims;
                    keys=nirs.getStimNames(data(i));
                    for k=1:length(keys)
                        if(~isa(data(i).stimulus(keys{k}),'nirs.design.StimulusVector'));
                            job.listOfStims{end+1,1}=keys{k};
                        end
                    end
                    data(i)=job.run(data(i));
                end
                
            end
            
        end
        
    end
    
end

