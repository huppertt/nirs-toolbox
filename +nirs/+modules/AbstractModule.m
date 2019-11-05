classdef AbstractModule
%% AbstractModule - Abstract class at the core of the pipeline design.
    
    properties
        name    = ''; % name of the module for convenience
        prevJob = []; % the module preceding this module in the pipeline
    end
    properties(Hidden=true);
        citation ='';
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
        
        function varargout = cite(obj)
            j=nirs.modules.pipelineToList(obj);
            cnt=1;
            for i=1:length(j)
                
                s=j{i}.citation;
                if(~iscellstr(s)); s=cellstr(s); end;
                for l=1:length(s)
                    ss='';
                    for k=1:100:length(s{1});
                        ss=sprintf('%s\n%s',ss,s{l}(k:min(k+99,length(s{l}))));
                    end;
                    if(length(ss)>0)
                        ss(1)=[];
                        
                    end
                    out(cnt)=cellstr(ss);
                    name{cnt}=j{i}.name;
                    cnt=cnt+1;
                end
            end
            [out,list]=unique(out);
            if(nargout==0)
                disp('Citations:');
                for i=length(out):-1:1
                    
                    if(~isempty(out{i}))
                        disp('-----------------')
                        disp(name{list(i)});
                        disp(out{i});
                    end
                end
            else
                varargout{1}=out;
            end
        end
        
        function out = run( obj, inputs )
            
            if(nargin<2)
                inputs=[];
            end
            
            % if no prev job execute and return result
            if isempty( obj.prevJob )
                if(iscell(inputs))
                    for i=1:length(inputs)
                          out{i} = obj.runThis( inputs{i} );
                    end
                else
               out = obj.runThis( inputs );
                end
                
            % else execute prev job first
            else
                 if(iscell(inputs))
                    for i=1:length(inputs)
                           out{i} = obj.runThis( obj.prevJob.run(inputs{i} ) );
                    end
                else
               out = obj.runThis( obj.prevJob.run(inputs ) );
                end
               
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

