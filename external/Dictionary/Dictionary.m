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
    %     toStruct - converts to a Matlab struct
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
        function tmp =toStruct(obj)
            tmp=struct;
            if(obj.count>0)
               for i=1:obj.count
                    tmp=setfield(tmp,obj.keys{i},obj.values{i});
                end
            end
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
                
                disp(['     ' obj.keys{i} '  :  ' str(:)']);
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
                    warning("%s does not exist in Dictionary",keys{k});
                end
            end
            
            n = length(out);
            if n == 1
                out = out{1};
            elseif n == 0
                out = [];
            end
        end

        function n = numArgumentsFromSubscript(obj,s,indexingContext)
           if indexingContext == matlab.mixin.util.IndexingContext.Expression
              n = 1;
           else
              if(iscell(s(1).subs))
                n = length(s(1).subs);
              else
                n = 1;
              end
           end
        end
        
        % assignment, i.e. dict('hello') = 1234
        function obj = subsasgn(obj,s,b)
            
            numSubRef=length(s);
            if(numSubRef>0)
                if strcmp(s(1).type,'.')
                    key = s.subs;
                    if isprop(obj,key)||ismethod(obj,key)
                        obj = builtin('subsasgn',obj,s);
                    else
                        s(1).type='()';
                        s(1).subs={s(1).subs};
                    end
                end
                
                if strcmp(s(1).type,'()')
                    % Assignment of dictionary item
                    for k=1:length(s(1).subs)
                        % for each item, assign it to new item b
                        newKey      = s(1).subs{k};
                        newValue    = b;
    
                        if(isempty(newValue)&&iscell(newValue))
                            % If assignment is empty field is removed from dictionary
                            % (if its a cell)
                            % allows easy deletion and empty strings
                            if(isempty(obj.get(newKey)))
                                obj.remove(newKey);
                            end
                        else
                            if(numSubRef>1) % If looking into subfields
                                dictItemToUpdate = obj.get( newKey ); %get the item from the dictionary to update
                                
                                if(isempty(dictItemToUpdate))
                                    error('Item "%s" does not exist in dictionary',newKey);
                                end
                                
                                dictItemToUpdate = subsasgn(dictItemToUpdate,s(2:end),b);
    
                                newValue=dictItemToUpdate;
                            end
                            obj = obj.put( newKey, newValue );
                        end
    
                    end
                else
               	    obj = builtin('subsasgn',obj,s,b);
                end
            end
        end
        
        % retrieval; i.e. dict('hello') returns 1234
        function [varargout] = subsref(obj,s)
            varargout=cell(1,1);
            out=cell(1,1);

            numSubRef=length(s);
            if(numSubRef>=1)
                switch s(1).type
                    case '()'
                        key = s(1).subs;
                        out{:} = obj.get( key );
                        % Return 1 output per key

                        if(numSubRef>1)
                            if(length(key)~=1)
                                out=out{1};
                            end
                            for j=1:length(key)
                                if(isempty(out{j}))
                                    error('Dictionary Item is empty!');
                                end
                                out{j} = builtin('subsref',out{j},s(2:end));
                            end
                        else
                            if(~iscell(out))
                                 out={out};
                            end
                        end
                    case '.'
                        key = s.subs;
                        if isprop(obj,key)||ismethod(obj,key)
                            out{:} = builtin('subsref',obj,s);
                        else
                            out{:} = obj.get( key );
                            if(isempty(out))
                                error('%s is not a property, method, or item of Dictionary object',key);
                            end
                            if(numSubRef>1)
                                key={key};
                                if(length(key)~=1)
                                    out=out{1};
                                end
                                for j=1:length(key)
                                    if(isempty(out{j}))
                                        error('Dictionary Item is empty!');
                                    end
                                    out{j} = builtin('subsref',out{j},s(2:end));
                                end
                            else
                                if(~iscell(out))
                                     out={out};
                                end
                            end
                        end
                    otherwise
                        out{:} = builtin('subsref',obj,s);
                end
            else
                out{:} = builtin('subsref',obj,s);
            end


            if(nargout==0)
                varargout=out;
            else
                for k=1:max(nargout,length(out)) % At least assign the first output
                    varargout{k}=out{k};
                end
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
                out = (length(keys) == length(unique(keys)));
        end
    end
    
    methods ( Access = private )
        % find index
        function [i, keyexists] = getindex( obj, key )
            if(iscell(key))
                key=key{1};
            end
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

