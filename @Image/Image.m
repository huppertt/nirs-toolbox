classdef Image
    %IMAGE This class holds segmented image volumes for forward and inverse
    % models.
    
    properties
        vol;        % image volume
        dim;        % spatial dimensions in mm
        origin;     % idx of (0,0,0)
        description; 
    end
    
    methods
        %% Constructor
        function obj = Image( vol, dim, origin )
            if nargin > 0, 
                obj.vol = vol;
                obj.dim = ones(1,ndims(vol));
                obj.origin = ones(1,ndims(vol));
            end
            if nargin > 1, obj.dim = dim; end 
            if nargin > 2, obj.origin = origin; end 
        end
        
        function obj = set.vol( obj, vol )
            assert( isnumeric(vol) || islogical(vol) )
            obj.vol = vol;
        end
        
        function obj = set.description( obj, description )
            assert( ischar( description ) || isempty( description ) )
          	obj.description = description;
        end
        
        function out = size( obj )
            out = size(obj.vol);
        end
        
    end
    
end

