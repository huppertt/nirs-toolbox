classdef AbstractModule
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name    = '';
        prevJob = [];
    end
    
    methods( Abstract )
       output   = execute       ( obj, input );
       options  = getOptions    ( obj );
       obj      = putOptions    ( obj, options );
    end
    
    methods
        function out = run( obj, input )
            % if no prev job execute and return result
            if isempty( obj.prevJob )
                out = obj.execute( input );
                
            % else execute prev job first
            else
                out = obj.execute( obj.prevJob.run( input ) );
            end
            
        end
        
    end
    
end

