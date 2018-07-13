classdef identity
% identity basis function for inverse models
    
    properties
        nVox;
    end
    properties( Dependent = true )
        fwd;
        inv;
    end
    methods
        function obj = identity(nVox)
            if nargin > 0;
                if(isa(nVox,'nirs.core.Mesh'))
                    nVox=size(nVox.nodes,1);
                end
                obj.nVox = nVox; 
            end
        end
        
        function out = get.fwd( obj)
            out=speye(obj.nVox,obj.nVox);
            
        end
         function out = get.inv( obj)
            out=speye(obj.nVox,obj.nVox);
            
        end
    end
    
end

