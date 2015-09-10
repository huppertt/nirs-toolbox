classdef gaussian
% identity basis function for inverse models
    
    properties
        Mtx;
    end
    properties( Dependent = true )
        fwd;
        inv;
    end
    methods
        function obj = gaussian(mesh,sigma,neighbors)
            if nargin < 1;
                error('Basis must be initialized with a mesh'); 
            end;
            
            if(nargin<2)
                sigma=10;
            end
            if(nargin<3)
                neighbors=5;
            end
            
            n=size(mesh.nodes,1);
            Mtx = zeros(n,n);
            for idx=1:size(mesh.faces,1)
                Mtx(mesh.faces(idx,1),mesh.faces(idx,2))=1;
                Mtx(mesh.faces(idx,1),mesh.faces(idx,3))=1;
                Mtx(mesh.faces(idx,2),mesh.faces(idx,3))=1;
                
                Mtx(mesh.faces(idx,2),mesh.faces(idx,1))=1;
                Mtx(mesh.faces(idx,3),mesh.faces(idx,1))=1;
                Mtx(mesh.faces(idx,3),mesh.faces(idx,2))=1;
            end
            for idx=2:neighbors
                Mtx =(Mtx*Mtx>0)*1;
            end
           
            X = repmat(mesh.nodes(:,1),1,n);
            Y = repmat(mesh.nodes(:,2),1,n);
            Z = repmat(mesh.nodes(:,3),1,n);
            
            X=X-X';
            Y=Y-Y';
            Z=Z-Z';
            Dist = X.^2 + Y.^2 + Z.^2;
            obj.Mtx = Mtx.*exp(-Dist/sigma^2);
            
        end
        
        function out = get.fwd( obj)
            out=obj.Mtx;
            
        end
         function out = get.inv( obj)
            out=obj.Mtx';
            
        end
    end
    
end

