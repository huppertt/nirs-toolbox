classdef AbstractModule
    %UNTITLED4 Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        name    = '';
        prevJob = [];
    end
    
    methods( Abstract )
       output   = runThis( obj, input );
    end
    
    methods
        function out = run( obj, input )
            % if no prev job execute and return result
            if isempty( obj.prevJob )
                out = obj.runThis( input );
                
            % else execute prev job first
            else
                out = obj.runThis( obj.prevJob.run( input ) );
            end
        end
        
        % option interface
        function out = options( obj, opts )
            if nargin == 1
                out = obj.getoptions();
            else
                out = obj.putoptions( opts );
            end
        end
    end
        
        
    methods( Access = private )
        function out = getoptions( obj )
            % we can inspect properties attributes with this
            mc = metaclass( obj );
            
            out = {};
            for i = 1:length( mc.PropertyList )
                p = mc.PropertyList(i);
                if strcmp(p.SetAccess, 'public') ...        % must be public
                        && ~p.Dependent  ...                % must not be dependent
                        && ~strcmp( p.Name, 'name' ) ...    % must not be "name"
                        && ~strcmp( p.Name, 'prevJob' ) ...     % must not be "prevJob"
                        && ~p.Hidden
                    
                    out = [out; {p.Name}];
                end
            end
            
         
        end
        
        function obj = putoptions( obj, opts )
            assert( isa( opts, 'Dictionary' ) )
            
            for i = 1:length(opts.keys)
                obj.(opts.keys{i}) = opts.values{i};
            end
        end
    end
        
    
end

