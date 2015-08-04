classdef ExportData < nirs.modules.AbstractModule
    %This function imports a variable from the workspace
    %   Detailed explanation goes here
    
    properties
        Output='raw';
    end
    
    methods

        function obj = ExportData( prevJob )
           obj.name = 'Export Data';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            assignin('base', obj.Output,data);
        end
    end
    
end