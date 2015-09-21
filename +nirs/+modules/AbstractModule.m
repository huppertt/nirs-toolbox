classdef AbstractModule
%% AbstractModule - Abstract class at the core of the pipeline design.
    
    properties
        name    = ''; % name of the module for convenience
        prevJob = []; % the module preceding this module in the pipeline
    end
    
    methods( Abstract )
        % every module must specify a runThis function which performs the
        % core functionality of that module
        output = runThis( obj, input );
    end
    
    methods
        % Every module automatically inherits the run function. This function
        % calls the previous module before performing its own
        % functionality.
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
            %% options - returns list of options or put options into module
            if nargin == 1
                out = obj.getoptions();
            else
                out = obj.putoptions( opts );
            end
        end
        
        function prop = javaoptions(obj)
            values = obj.getoptions();
            for idx=1:length(values)
                prop(idx)=javatypes(class(obj.(values{idx})));
                set(prop(idx),'Name',values{idx},'Value',obj.(values{idx}));
                set(prop(idx),'Category','Misc');
                set(prop(idx),'Description','');
            end
            if(~exist('prop')), prop=[]; end;
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

