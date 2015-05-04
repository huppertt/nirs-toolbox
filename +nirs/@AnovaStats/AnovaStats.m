classdef AnovaStats
    %ANOVASTATS Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        names	% variable names
        F       % F-stats
        df1     % degrees of freedom 1
        df2     % degrees of freedom 2
        probe   % probe geometry
        
        pcrit = 0.05;
    end
    
    properties ( Dependent = true )
        Fcrit
        p
    end
    
    methods
        function draw( obj )
           pass 
        end
    end
    
end