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
        function obj = gaussian(mesh,reducer)
            if nargin < 2;
                error('Basis must be initialized with a mesh'); 
            end;
            
            f=figure('visible','off');
            pv=patch('Faces',mesh.faces,'Vertices',mesh.nodes);
            pv=reducepatch(pv,reducer);
            close(f);
            nV=size(mesh.nodes,1);
            Mtx=zeros(nV,size(pv.vertices,1));
            for idx=1:size(Mtx,2)
                Mtx(:,idx)=(sum((mesh.nodes-ones(nV,1)*pv.vertices(idx,:)).^2,2));
            end
            sigma2=median(min(Mtx,[],2));
            
            Mtx=exp(-Mtx/sigma2);
            Mtx=Mtx./(max(Mtx,[],2)*ones(1,size(Mtx,2)));
            Mtx=Mtx.*(Mtx>1E-3);
           % Mtx=orth(Mtx);
            obj.Mtx=Mtx.*(abs(Mtx)>1E-3);
        end
        
        function out = get.fwd( obj)
            out=obj.Mtx;
            
        end
         function out = get.inv( obj)
            out=obj.Mtx';
            
        end
    end
    
end

