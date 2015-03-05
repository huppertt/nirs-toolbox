classdef RemoveStimless < nirs.functional.AbstractModule
    %UNTITLED6 Summary of this class goes here
    %   Detailed explanation goes here
    
    
    methods
        function obj = RemoveStimless( prevJob )
           obj.name = 'Remove Files w/o Stim';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = execute( obj, data )
            for i = 1:length(data)
                if isempty(data(i).stimulus.keys)
                    lst(i) = false;
                else
                    lst(i) = true;
                end
            end
            
            data = data(lst);
            
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
    end
    
end

