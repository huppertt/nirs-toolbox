classdef HashTable
    % This is not a real hash table, but will act like one with reasonable
    % performance for small tables. Matlab's builtin containers.Map is a
    % handle class, which causes huge problems when used as a property in
    % a value class.  This will work until I can find a suitable
    % replacement.
    
    properties (SetAccess = private)
        keys = {};
        values = {};
    end
    
    properties( Dependent = true )
        count
    end
    
    methods
        
        function obj = HashTable( keys, values )
            if nargin == 2
                assert( iscellstr(keys) )
                obj.keys = keys;
                obj.values = values;
            elseif nargin == 1
                error('Constructor takes zero or two arguments.')
            end
        end
            
        function count = get.count( obj )
            count = length( obj.keys );
        end
        
        function out = subsref(obj,s)
            if length(s) == 1 && strcmp(s.type,'()')
                if ischar(s.subs{1})
                    s = s.subs;
                elseif iscell( s.subs{1} )
                    assert(length(s.subs) == 1)
                    s = s.subs{1};
                end
                
                for i = 1:length(s)
                    lst = strcmp(obj.keys, s{i});
                    if any(lst)
                        out{i,1} = obj.values{lst};
                    else
                        out{i,1} = [];
                    end
                end

                if length(out) == 1
                    out = out{1};
                end
            else
                out = builtin('subsref',obj,s);
            end
        end
        
        function obj = subsasgn(obj,s,b)
            if strcmp(s.type,'()')
                if ischar(s.subs{1})
                    s = s.subs;
                elseif iscell( s.subs{1} )
                    assert(length(s.subs) == 1)
                    s = s.subs{1};
                end
                
                if ~iscell(b)
                    b = {b};
                end
                
                assert( length(s)==length(b) )
            
                for i = 1:length( s )
                    lst = strcmp(obj.keys, s{i});
                    if any(lst)
                        obj.values{lst} = b{i};
                    else
                        obj.keys{end+1} = s{i};
                        obj.values{end+1} = b{i};
                    end
                end
                
            else
               	obj = builtin('subsasgn',obj,s,b);
            end
        end
        
        function obj = delete( obj, keys )
            if ischar(keys)
                keys = {keys};
            end
            
            for i = 1:length( keys )
                lst = strcmp( keys{i}, obj.keys );
                obj.keys = obj.keys(~lst);
                obj.values = obj.values(~lst);
            end
        end
        
        function out = iskey( obj, keys )
            if ischar( keys )
               keys = {keys}; 
            end
            
            for i = 1:length( keys )
                lst = strcmp( keys{i}, obj.keys );
                if any(lst)
                    out(i) = true;
                else
                    out(i) = false;
                end
            end
        end
        
        function out = isempty( obj )
           if obj.count == 0
               out = true;
           else
               out = false;
           end
        end
        
    end
    
end

