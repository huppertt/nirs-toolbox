classdef Image
    %IMAGE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        dim;
        origin;
        refPts;
        volume;
    end
    
    methods
        function obj = Image( newVolume, varargin )
            obj.volume = newVolume;
            obj.origin = ones( ndims(obj.volume),1 ); %default origin [1 1 1]
            obj.dim = ones( ndims(obj.volume),1 );   %default dimensions [1 1 1] mm
            
            if nargin > 1
                obj.origin = varargin{2};
            end
            if nargin > 2
                obj.dim = varargin{3};
            end
            if nargin > 3
                obj.refPts = varargin{4};
            end
            if nargin > 4
                error('Too many input arguments for Image constructor.')
            end
        end
        
%         function obj = set.dim( obj, newDim )
%             if (length(newDim) == 3 ||length(newDim) == 4)
%                 if iscolumn(newDim)
%                     obj.dim = newDim;
%                 elseif isrow(newDim)
%                     obj.dim = newDim';
%                 end
%             else
%                 error('Image dimensions must be 3x1 or 4x1.')
%             end
%         end
        
        function obj = set.volume( obj, newVolume )
            if isnumeric( newVolume ) || islogical( newVolume )
                obj.volume = newVolume;
            else
                error( 'Volume must be numeric or logical.' )
            end
        end
        
        function dispImage( obj )
            view_nii( make_nii( obj.volume ) );
        end

%         function out = subsref(obj,s)
%             if strcmp(s(1).type,'()')
%                 out = obj.volume(s(1).subs{:});
%             else
%                 try
%                     out = builtin('subsref',obj,s);
%                 catch
%                     builtin('subsref',obj,s);
%                 end
%             end
%         end
%                 
%         function obj = subsasgn(obj,s,b)
%             if strcmp(s.type,'()')
%                 obj.volume(s.subs{:}) = b;
%             else
%                 obj = builtin('subsasgn',obj,s,b);
%             end
%         end 
    end
    
end

