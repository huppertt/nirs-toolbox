classdef ImportData < nirs.modules.AbstractModule
    %This function imports a variable from the workspace
    %   Detailed explanation goes here
    
    properties
        Input='raw';
    end
    
    methods

        function obj = ImportData( prevJob )
           obj.name = 'Import Data';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            data = evalin('base', obj.Input);
        end
    end
    
end