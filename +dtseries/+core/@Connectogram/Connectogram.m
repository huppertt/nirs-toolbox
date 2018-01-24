classdef Connectogram
    properties
        DataMatrix
        LabelsTable
    end
    
    methods
        function obj = Connectogram(data,labels)
            if(nargin<1)
                DataMatrix=[];
            end
            if(nargin<2)
                LabelsTable=table();
            end
        end
    
        function h = draw(obj)
            f=figure;
            
            
            
        end
        
    end
    
    
    
end