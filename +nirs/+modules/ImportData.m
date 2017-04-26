classdef ImportData < nirs.modules.AbstractModule
    %This function imports a variable from the workspace
    %   Detailed explanation goes here
    
    properties
        Input='raw';
        override=false;
    end
    
    methods

        function obj = ImportData( prevJob )
           obj.name = 'Import Data';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            if(isempty(data) | obj.override)
                data = evalin('base', obj.Input);
            else
                disp('Import data module skipped: Data provided');
            end
        end
        function prop = javaoptions(obj)
           prop = javaoptions@nirs.modules.AbstractModule(obj);
           set(prop(1),'Category','ImportData',...
               'Description','Name of the matlab variable to import');
           
        end
        
        
    end
    
end