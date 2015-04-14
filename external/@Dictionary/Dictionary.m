% Copyright (c) 2015, Jeffrey W Barker (jwb52@pitt.edu)
% All rights reserved.
% 
% Redistribution and use in source and binary forms, with or without
% modification, are permitted provided that the following conditions
% are met:
% 
% 1. Redistributions of source code must retain the above copyright
% notice, this list of conditions and the following disclaimer.
% 
% 2. Redistributions in binary form must reproduce the above copyright
% notice, this list of conditions and the following disclaimer in the
% documentation and/or other materials provided with the distribution.
% 
% THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
% "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
% LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
% A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT
% HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT,
% INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING,
% BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS
% OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
% AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
% LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY
% WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
% POSSIBILITY OF SUCH DAMAGE.

classdef Dictionary
    
    properties (SetAccess = private)
        keys
        values
    end
    
    properties( Dependent = true )
        count
    end
    
    properties ( Access = private )
        indices;
        TABLE_SIZE = uint32(1024);
    end
    
    methods
        
        % constructor
        function obj = Dictionary( keys, vals )
            
            obj.keys    = {};
            obj.values  = {};
            
            if nargin == 2
                assert( iscellstr(keys) ....
                    && iscell(vals) ....
                    && length(keys) == length(vals) ....
                    )
                
                obj.TABLE_SIZE  = 2 * length(keys);
                obj.keys        = keys;
                obj.values      = vals;
            elseif nargin == 1
                error('Constructor takes zero or two arguments.')
            end
            
            obj = obj.rehash();
        end
        
        % update with list of keys and vals
        function obj = update(obj, keys, vals)
            assert( iscellstr(keys) && length(keys)==length(vals) )                
            for i = 1:length(keys)
                obj.put(keys{i},vals{i});
            end
        end
        
        % number of items in dictionary
        function count = get.count( obj )
            count = length(obj.keys);
        end
        
        % delete items
        function obj = delete( obj, keys )
            if ischar(keys)
                keys = {keys};
            end
            
            for k = 1:length(keys)
               [i, keyexists] = obj.getindex(keys{k});
               if keyexists
                   idx = obj.indices(i);
                   obj.keys(idx) = [];
                   obj.values(idx) = [];

                   lst = obj.indices > idx;
                   obj.indices(lst) = obj.indices(lst) - 1;
               end
            end
            
        end
        
        % check if keys exists
        function out = iskey( obj, key )
            [~,keyexists] = obj.getindex(key);
            out = keyexists;
        end
        
        % check if empty
        function out = isempty( obj )
           out = obj.count == 0;
        end
        
        % assignment, i.e. dict('hello') = 1234
        function obj = subsasgn(obj,s,b)
            if strcmp(s.type,'()')
                assert( ischar(s.subs{1}) )
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
                assert( ischar(s.subs{1}) )
                key = s.subs{1};
                out = obj.get( key );
            else
                out = builtin('subsref',obj,s);
            end
        end
        
        function obj = resize( obj, N )
            assert( N < uint32(2^32-1) )
   
            % resize table
            obj.TABLE_SIZE   = N;
            
            % rehash indices
            obj = obj.rehash();
        end
    end
    
    methods ( Access = private )
        function h = hash( obj, s )
            % this is faster than anything that can be 
            % implemented in pure matlab code
            h = typecast(java.lang.String(s).hashCode(),'uint32');
            h = mod(h(2), obj.TABLE_SIZE) + 1;
        end
        
        % insert new items
        function obj = put( obj, newKey, newValue )
            if (obj.count + 1) > idivide(obj.TABLE_SIZE, uint32(2))
               obj = obj.resize( 2 * obj.TABLE_SIZE ); 
            end

            [i, keyexists] = obj.getindex( newKey );
            
            if keyexists % key already exists
                idx = obj.indices(i);
                obj.values{idx} = newValue;
            else
                % add index
                obj.indices(i) = length(obj.keys)+1;
                
                % append keys and values
                obj.keys    {end+1} = newKey;
                obj.values  {end+1} = newValue;
            end
        end
        
        % get items
        function out = get( obj, key )
            [i, keyexists] = obj.getindex(key);
            
            if keyexists
                idx = obj.indices(i);
                out = obj.values{idx};
            else
                out = [];
            end
        end
        
        % find index
        function [i, keyexists] = getindex( obj, key )
            
            i = obj.hash( key );
            while obj.indices(i) > 0 ...                    % entry full
                && ~strcmp(obj.keys{obj.indices(i)}, key)  	% key doesn't match
             	i = mod(i, obj.TABLE_SIZE) + 1;
            end
            
            keyexists = obj.indices(i) > 0;
        end
        
        % rehash indices
        function obj = rehash( obj )
            obj.indices = zeros(obj.TABLE_SIZE,1,'uint32');
            for k = 1:length(obj.keys)
                   i = obj.getindex(obj.keys{k});
                   obj.indices(i) = k;
            end
        end
    end
    
end

