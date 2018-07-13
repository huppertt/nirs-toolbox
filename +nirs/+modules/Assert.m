classdef Assert < nirs.modules.AbstractModule
%% Throws a gentle error if condition is not met
%
% Options:
%     condition string or function handle evaluated
    
    properties
        condition=@(data)isa(data.probe,'nirs.core.Probe1020');
        throwerror=false;
    end
    
    methods
        function obj = Assert( prevJob )
           obj.name = 'Assert for job errors';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            if(~iscell(obj.condition))
                obj.condition={obj.condition};
            end
            files={};
            errormsg={};
            for i = 1:numel(data)
                files{i}=data(i).description;
                for j=1:length(obj.condition)
                    
                    try;
                        bool(i,j)=feval(obj.condition{j},data(i));
                    catch
                        bool(i,j)=false;
                        errormsg{i}=lasterr;
                    end
                end
                
            end
            if(~all(bool(:)) )
                for j=1:length(obj.condition)
                    disp(func2str(obj.condition{j}))
                    for i = 1:numel(data)
                        if(bool(i,j))
                            disp(['    ' num2str(i) ':' data(i).description ' PASSED']);
                        else
                            disp(['    ' num2str(i) ':' data(i).description ' FAILED']);
                        end
                    end
                end
                
                if(obj.throwerror)
                    error('nirs.modules.Assert failed');
                end
            end
        end
    end
    
end

