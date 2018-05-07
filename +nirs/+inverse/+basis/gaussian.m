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
                neighbors=3;
            end
            
            n=size(mesh.nodes,1);
            if(neighbors>0)
                
                if(isempty(mesh.elems))
                    i=[mesh.faces(:,1); mesh.faces(:,1); mesh.faces(:,2); mesh.faces(:,2);...
                        mesh.faces(:,3); mesh.faces(:,3)];
                    j=[mesh.faces(:,2); mesh.faces(:,3); mesh.faces(:,3); mesh.faces(:,1);...
                        mesh.faces(:,1); mesh.faces(:,2)];
                else
                    i=[mesh.elems(:,1); mesh.elems(:,1); mesh.elems(:,1); mesh.elems(:,2); mesh.elems(:,2); mesh.elems(:,2);...
                        mesh.elems(:,3); mesh.elems(:,3); mesh.elems(:,3); mesh.elems(:,4); mesh.elems(:,4); mesh.elems(:,4)];
                    j=[mesh.elems(:,2); mesh.elems(:,3); mesh.elems(:,4); mesh.elems(:,3); mesh.elems(:,1); mesh.elems(:,4);...
                        mesh.elems(:,1); mesh.elems(:,2); mesh.elems(:,4); mesh.elems(:,1); mesh.elems(:,2); mesh.elems(:,3)];
                    
                end
                Mtx = sparse([i' 1:n],[j' 1:n],ones(size(i,1)+n,1),n,n);
%                 
%                 Mtx(sub2ind([n n],i,j))=1;
%                 Mtx(sub2ind([n n],1:n,1:n))=1;
%                 Mtx=sparse(Mtx);
                %
                %             for idx=1:size(mesh.faces,1)
                %                 Mtx(mesh.faces(idx,1),mesh.faces(idx,2))=1;
                %                 Mtx(mesh.faces(idx,1),mesh.faces(idx,3))=1;
                %                 Mtx(mesh.faces(idx,2),mesh.faces(idx,3))=1;
                %
                %                 Mtx(mesh.faces(idx,2),mesh.faces(idx,1))=1;
                %                 Mtx(mesh.faces(idx,3),mesh.faces(idx,1))=1;
                %                 Mtx(mesh.faces(idx,3),mesh.faces(idx,2))=1;
                %             end
                for idx=2:neighbors
                    Mtx =(Mtx*Mtx>0)*1;
                end
            end
            
            X = repmat(mesh.nodes(:,1),1,n);
            Y = repmat(mesh.nodes(:,2),1,n);
            Z = repmat(mesh.nodes(:,3),1,n);
            
            X=X-X';
            Y=Y-Y';
            Z=Z-Z';
            Dist = X.^2 + Y.^2 + Z.^2;
            if(neighbors>0)
                obj.Mtx = Mtx.*exp(-Dist/sigma^2);
            else
                obj.Mtx = exp(-Dist/sigma^2);
                obj.Mtx = obj.Mtx.*(obj.Mtx>exp(-3));
            end
            
        end
        
        function out = get.fwd( obj)
            out=obj.Mtx;
            
        end
        function out = get.inv( obj)
            out=obj.Mtx';
            
        end
    end
    
end

