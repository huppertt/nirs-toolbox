classdef RunMatlabCode < nirs.modules.AbstractModule
    %%  RunMatlabCode - This is a generic wrapper function to run an feval statement
    % on a data file as part of the pipeline
    %
    % Options:
    %     FunctionHandle - a filename or function handle to call
    %                    must be of the form Data = fcn(Data)
    
    properties
        FunctionHandle   % filename or function handle to call
        rungroup;
    end
    
    methods
        
        function obj =  RunMatlabCode( prevJob )
            obj.name = 'Pass data into generic matlab function';
            obj.rungroup=false;
            if nargin > 0
                obj.prevJob = prevJob;
            end
        end
        
        function data = runThis( obj, data )
            % columns of the table that arent varToMatch
            if(nargin(obj.FunctionHandle)==0)
                feval(obj.FunctionHandle);
            else
                if(~obj.rungroup)
                    for i = 1:numel(data)
                        data(i)=feval(obj.FunctionHandle,data(i));
                        if(isa(data,'nirs.core.Data') && isa(data(i).data,'nirs.core.ChannelStats'))
                            data=data.data;
                        end
                    end
                else
                    data=feval(obj.FunctionHandle,data);
                end
            end
        end
        
    end
end
