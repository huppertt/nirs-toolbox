classdef freesurfer_wavelet
% identity basis function for inverse models
    
    properties
        ntypes=2;
        nlevels;
        vert_per_level;
        mesh_verts;
        WlSynthesisMtx;
    end
    properties( Dependent = true )
        fwd;
        inv;
    end
    methods
        function obj = freesurfer_wavelet(nLevels,nLevelsmesh,ntypes)
            
             p = fileparts( which('nirs.modules.ImageReconMFX') );
             W = load([p filesep 'wavelet_matrix.mat']);
             if(nargin<1)
                 nLevels=length(W.template.vertex);
             else
                 nLevels=min(nLevels,length(W.template.vertex));
             end
            if(nargin>2)
                obj.ntypes=ntypes;
            end
            if(nargin<2) 
                nLevelsmesh=6;
            end
            
             for i=1:nLevels
                 obj.vert_per_level(i)=length(W.template.vertex{i});
             end
             obj.mesh_verts=length(W.template.vertex{nLevelsmesh});
            
             if nargin > 0;
                obj.nlevels = length(obj.vert_per_level); 
             end
            obj.WlSynthesisMtx = W.WlSynthesisMtx; 
            
        end
        
        function out = get.fwd( obj)
            out=obj.WlSynthesisMtx(1:obj.mesh_verts,1:obj.vert_per_level(obj.nlevels));
            for i=2:obj.ntypes
                out = blkdiag(out,obj.WlSynthesisMtx(1:obj.mesh_verts,1:obj.vert_per_level(obj.nlevels)));
            end
        end
         function out = get.inv( obj)
                out=obj.WlSynthesisMtx(1:obj.mesh_verts,1:obj.vert_per_level(obj.nlevels)); 
                for i=2:obj.ntypes
                    out = blkdiag(out,obj.WlSynthesisMtx(1:obj.mesh_verts,1:obj.vert_per_level(obj.nlevels)));
                end
        end
    end
    
end

