classdef Mesh

    properties
        nodes
        faces
        elems
    end
    
    methods
        function obj = Mesh( n, f, e )
           if nargin > 0, obj.nodes = n; end
           if nargin > 1, obj.faces = f; end
           if nargin > 2, obj.elems = e; end
        end
        
        function h = draw( obj, values, vmax, thresh, cmap )
            if nargin < 5, cmap     = []; end  
            if nargin < 4, thresh   = []; end
            if nargin < 3, vmax     = []; end
            if nargin < 2, values   = []; end
            
            h = nirs.plotmesh( obj.nodes, obj.faces, ...
                	values, vmax, thresh, cmap );
        end
    end
    
end

