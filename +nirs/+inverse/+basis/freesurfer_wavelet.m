classdef freesurfer_wavelet
% identity basis function for inverse models
    
    properties
        nlevels;
        vert_per_level;
        WlSynthesisMtx;
    end
    properties( Dependent = true )
        fwd;
        inv;
    end
    methods
        function obj = freesurfer_wavelet(nLevels)
            
%             
             p = fileparts( which('nirs.modules.ImageReconMFX') );
             W = load([p filesep 'wavelet_matrix.mat']);
             
             for i=1:length(W.template.vertex)
                 obj.vert_per_level(i)=length(W.template.vertex{i});
             end
      
             if nargin > 0;
                obj.nlevels = length(obj.vert_per_level); 
             end
            obj.WlSynthesisMtx = W.WlSynthesisMtx; 
            
   
    
    
            
        end
        
        function out = get.fwd( obj)
            out=obj.WlSynthesisMtx(:,1:obj.vert_per_level(obj.nlevels));
            out = blkdiag(out,out);
                     
        end
         function out = get.inv( obj)
                out=obj.WlSynthesisMtx(:,1:obj.vert_per_level(obj.nlevels));  
                out = blkdiag(out,out);
        end
    end
    
end

