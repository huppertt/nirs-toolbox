classdef AbstractJob
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name = '';
        prevJob = [];
    end
    
    methods( Abstract )
       output   = execute       ( obj, input );
       options  = getOptions    ( obj );
       obj      = putOptions    ( obj, options );
    end
    
    methods
        function out = run( obj, input )
            if isempty( obj.prevJob )
                out = obj.execute( input );
            else
                out = obj.execute( obj.prevJob.run( input ) );
            end
            
        end
        
    end
    
end

