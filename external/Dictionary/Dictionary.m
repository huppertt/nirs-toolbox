classdef Dictionary
    %% DICTIONARY - This is an implementation of a Hash Table/Dictionary that
    %               is a value class (can be copied with = operator).
    % 
    % Properties: 
    %     keys   - cell array of keys (these can be any data type)
    %     values - cell array of values for each key
    %     count  - number of items in the dictionary
    % 
    % Methods:
    %     iskey   - checks if key exists
    %     isempty - checks if Dictionary is empty
    %     remove  - removes an item from the Dictionary given a key
    %     hash    - (static) returns hash of a key
    % 
    % Example:
    %     >> d = Dictionary();
    %     >> d('a') = 1;
    %     >> d(@(x) x*x) = 2;
    %     >> d('a')
    % 
    %     ans =
    % 
    %          1
    % 
    %     >> d(@(x) x*x)
    % 
    %     ans =
    % 
    %          2
    
    properties (SetAccess = private)
        keys    = {} % cell array of keys (these can be any data type)
        values  = {} % cell array of values for each key
    end
    
    properties( Dependent = true )
        count % number of items in the dictionary
    end
    
    properties ( Access = private )
        indices;
        TABLE_SIZE = uint64(1024);
    end
    
    methods
        function obj = Dictionary( keys, vals )
            %% Dictionary - Creates a Dictionary object.
            % 
            % Args:
            %     keys - cell array of keys
            %     vals - cell array of vals
            %     
            % Examples:
            %     % empty dictionary
            %     d = Dictionary()
            %     
            %     % intitialized dictionary
            %     d = Dictionary({'a', 'b','c'}, {1, 2, 'hello'});
            
            % check java classpath
            
%             dictLoc = fileparts(which('Dictionary'));
%             clsPth  = javaclasspath('-dynamic');
%             if ~any( strcmp(dictLoc, clsPth) )
%                 javaaddpath( dictLoc );
%             end
            
            % if key/value pairs provide add to dict
            if nargin == 2
                assert( length(keys)==length(vals) ...
                    && iscell(vals) ...
                    && iscell(keys) ...
                    && Dictionary.areUniqueKeys(keys) )  

                obj.keys        = keys;
                obj.values      = vals;
            elseif nargin == 1
                error('Constructor takes zero or two arguments.')
            end
            
            obj = obj.rehash();
        end
        
        function count = get.count( obj )
            %% count - returns the number of items in dictionary
            count = length(obj.keys);
        end
        
        function obj = remove( obj, keys )
            %% remove - remove items from dictionary
            % 
            % Args:
            %     keys - either a single key or cell array of keys
            
            if ischar(keys)
                keys = {keys};
            end
            
            for k = 1:length(keys)
               [i, keyexists] = obj.getindex(keys{k});
               if keyexists
                   idx = obj.indices(i);
                   
                   obj.keys(idx)    = [];
                   obj.values(idx)  = [];
                   obj.indices(i)   = 0;

                   lst = obj.indices > idx;
                   obj.indices(lst) = obj.indices(lst) - 1;
               end
            end
            
        end
        
        function out = iskey( obj, key )
            %% iskey - returns true if key exists
            
            [~,keyexists] = obj.getindex(key);
            out = keyexists;
        end
        
        function disp(obj)
           % tt=table(obj.values{:},'VariableNames',obj.keys);
           if(obj.count>0)
            disp(['Dictionary Class Containing:'])
            for i=1:obj.count
                if(isnumeric(obj.values{i}))
                    str=num2str(obj.values{i});
                elseif(ischar(obj.values{i}))
                    str=obj.values{i};
                else
                    str=class(obj.values{i});
                end
                
                disp(['     ' obj.keys{i} '  :  ' str]);
            end
            %disp(tt);
           else
               disp(['Empty Dictionary Class'])
           end
            
            
        end
        
        
        function out = isempty( obj )
            %% isempty - returns true if dictionary is empty
            out = obj.count == 0;
        end
        
        function obj = resize( obj, N )
            %% resize - resize the  dictionary table
            %
            % Args:
            %     N - the size of the new table (should be 2-3 times
            %         the expected number of items to be put in the dictionary
            %
            % Notes: The Dictionary will resize itself as needed so
            %        manually resizing is not typically necessary
            
            assert( N < uint64(2^32) )
   
            % resize table
            obj.TABLE_SIZE   = N;
            
            % rehash indices
            obj = obj.rehash();
        end
    end
       
    methods (Hidden = true)
        % put new items
        function obj = put( obj, newKeys, newVals )
            
            if ~iscell(newKeys)
                newKeys = {newKeys};
            end
            
            if ~iscell(newVals)
                newVals = {newVals};
            end
            
            for k = 1:length( newKeys )
                [i, keyexists] = obj.getindex( newKeys{k} );

                if keyexists % key already exists
                    idx = obj.indices(i);
                    obj.values{idx} = newVals{k};
                else
                    % resize table if necessary
                    if (obj.count + 1) > obj.TABLE_SIZE / uint64(2)
                       obj = obj.resize( 3 * obj.TABLE_SIZE ); 
                    end
                    
                    % add index
                    obj.indices(i) = obj.count + 1;

                    % append keys and values
                    obj.keys    {end+1} = newKeys{k};
                    obj.values  {end+1} = newVals{k};
                end
            end
        end
        
        % get items
        function out = get( obj, keys )
            if ~iscell(keys)
                keys = {keys};
            end
            
            for k = 1:length( keys )
                [i, keyexists] = obj.getindex(keys{k});

                if keyexists
                    idx = obj.indices(i);
                    out{k} = obj.values{idx};
                else
                    out{k} = [];
                end
            end
            
            n = length(out);
            if n == 1
                out = out{1};
            elseif n == 0
                out = [];
            end
        end
        
        % assignment, i.e. dict('hello') = 1234
        function obj = subsasgn(obj,s,b)
            if strcmp(s.type,'()')
                % assert( ischar(s.subs{1}) )
                newKey      = s.subs{1};
                newValue    = b;

                obj = obj.put( newKey, newValue );
            else
               	obj = builtin('subsasgn',obj,s,b);
            end
        end
        
        % retrieval; i.e. dict('hello') returns 1234
        function out = subsref(obj,s)
            if length(s) == 1 && strcmp(s.type,'()')
                % assert( ischar(s.subs{1}) )
                key = s.subs{1};
                out = obj.get( key );
            else
                out = builtin('subsref',obj,s);
            end
        end
    end
    
    methods ( Static )
        function [h, b] = hash( key )
            %% hash - returns the hash of a key
             b = getByteStreamFromArray(key);
%             h = uint64(typecast(int32(MyHashLib.jenkinsHash(b)),'uint32'));
            
            % this can be used instead, but it is slower
            h = typecast(java.lang.String(b).hashCode(),'uint32');
            h = uint64(h(2));
        end
        
        function out = areUniqueKeys( keys )
            %% areUniqueKeys - returns true if the list of keys are unique
%             b = {};
%             for i = 1:length( keys )
%                b{i} = cast(getByteStreamFromArray(keys{i}),'uint32');
%             end
%             out = (length(b) == length(unique(b)));
                out = (length(keys) == length(unique(keys)));
        end
    end
    
    methods ( Access = private )
        % find index
        function [i, keyexists] = getindex( obj, key )
             b = getByteStreamFromArray(key);
%             i = uint64(typecast(int32(MyHashLib.jenkinsHash(b)),'uint32'));
%             
            i = typecast(java.lang.String(b).hashCode(),'uint32');
            i = uint64(i(2));
            
            i = mod(i, obj.TABLE_SIZE) + uint64(1);
            
            % while full and keys don't match
            while obj.indices(i) && ~isequal(b, getByteStreamFromArray(obj.keys{obj.indices(i)}))
                    % increment index; 
                    % wrap to beginning if necessary
                    if i == obj.TABLE_SIZE
                        i = uint64(1);
                    else
                        i = i + uint64(1);
                    end
            end
            
             keyexists = logical( obj.indices(i) );
        end
        
        % rehash indices
        function obj = rehash( obj )
            if(verLessThan('matlab', '9.2'))
                obj.TABLE_SIZE  = max(uint64(1024), uint64(4 * length(obj.keys)));
            else
                % The max(uint64) compltetly locks up matlab on vr > 2017
                % (Matlab 9.2)
                obj.TABLE_SIZE  = max(1024, 4 * length(obj.keys));
            end
            obj.indices = zeros(obj.TABLE_SIZE,1,'uint32');
            for k = 1:length(obj.keys)
                   i = obj.getindex(obj.keys{k});
                   obj.indices(i) = k;
            end
        end
    end
    
end

