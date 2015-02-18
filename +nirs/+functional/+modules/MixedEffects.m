classdef MixedEffects < nirs.functional.AbstractModule
    %UNTITLED3 Summary of this class goes here
    %   Detailed explanation goes here
  
    properties
        formula = 'beta ~ groupID';
    end
    
    methods

        function obj = MixedEffects( prevJob )
           obj.name = 'Mixed Effects Model';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function groupStats = execute( obj, subjStats )
            % assemble table
            for i = 1:length(subjStats)
                
                %% TODO
                
            end
            
            % call lme package
            
            % return stats
            
        end
        
        function options = getOptions( obj )
            options = [];
        end
           
        function obj = putOptions( obj, options )
        end
        
    end
    
end

