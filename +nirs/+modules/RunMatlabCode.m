classdef RunMatlabCode < nirs.modules.AbstractModule
%%  RunMatlabCode - This is a generic wrapper function to run an feval statement 
% on a data file as part of the pipeline 
% 
% Options: 
%     FunctionHandle - a filename or function handle to call
%                    must be of the form Data = fcn(Data)

    properties
        FunctionHandle   % filename or function handle to call

    end
    
    methods

        function obj =  RunMatlabCode( prevJob )
           obj.name = 'Pass data into generic matlab function';
           
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            % columns of the table that arent varToMatch
                
            for i = 1:length(data)
                data(i)=feval(obj.FunctionHandle,data(i));             
            end
        end
    
    end
end
