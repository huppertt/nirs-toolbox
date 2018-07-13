classdef ChannelStatsBootStrap
    properties
        pipeline
        samplefunction
        output = 'beta'
        onesided=false;
    end
    
    properties (SetAccess = protected)
        samples
    end
    
    methods
        function obj = ChannelStatsBootStrap(pipeline, samplefunction)
            if(nargin<1)
                obj.pipeline =nirs.modules.default_modules.single_subject;
                obj.samplefunction = @nirs.testing.simData;
            else
                obj.pipeline = pipeline;
            end
            if nargin > 1
                obj.simfunc = simfunc;
            end
        end
        
        function obj = run(obj,iter)
            for i=1:iter
                disp(['iter ' num2str(i) ' of ' num2str(iter)]);
                data = feval(obj.samplefunction);
                
                if(~isempty(strfind(class(obj.pipeline),'.modules')))
                    S = obj.pipeline.run(data);
                else
                    S =feval(obj.pipeline,data);
                end
                vals = vertcat(S.(obj.output));
                obj.samples=vertcat(obj.samples,vals(:));
                
            end
        
        end
        
        function p = pval(obj,values)
            p=ones(size(values));
            for i=1:length(values(:));
                if(obj.onesided)
                    p(i)=length(find(obj.samples(:)>=values(i)));
                else
                    if(values(i)<=0)
                        p(i)=length(find(abs(obj.samples(:))>=abs(values(i))));
                    end
                end
            end
            p=p/prod(size(obj.samples));
        end
        
    end
end