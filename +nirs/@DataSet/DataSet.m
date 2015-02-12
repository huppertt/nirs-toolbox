classdef DataSet
    %UNTITLED8 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        data
        stimulusConverters
    end
    
    properties( Dependent = true )
        conditionNames
        
        demographicVariables
    end
    
    methods
        function renameCondition( obj, oldName, newName )
            
        end
        
        function mergeConditions( obj, condNames, newName )
            
        end
        
        
    end
    
end

