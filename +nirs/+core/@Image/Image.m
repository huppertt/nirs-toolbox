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
        
        function mesh = convertFEMmesh(obj)
            %This converts a volume to a FEM mesh using iso2mesh
            try
                iso2meshver;
                disp('Using Iso2Mesh.  Please cite:');
                disp('Qianqian Fang and David Boas, "Tetrahedral mesh generation from volumetric binary')
                disp('and gray-scale images," Proceedings of IEEE International Symposium on Biomedical ');
                disp('Imaging 2009, pp. 1142-1145, 2009');
                
            catch
               
                disp('Please download the iso2mesh package from:');
                disp('http://iso2mesh.sourceforge.net');
                disp('and/or add to the matlab path');
                 error('Cannot find Iso2Mesh on Matlab Path');
                
            end;
            
            [n,e,f]=v2m(obj.vol,.01*max(obj.vol(:)),1,max(obj.dim));
            
            %Correct the scaling
            n=n.*(ones(size(n,1),1)*obj.dim);
            n=n+ones(size(n,1),1)*obj.origin;
            
            mesh = nirs.core.Mesh(n,f(:,1:end-1),e(:,1:end-1));
            mesh.regions=ones(size(n,1),1);
        end
        
        function mesh = convertBEMmesh(obj)
            %This converts a volume to a BEM (surface) mesh using iso2mesh
            try
                iso2meshver;
                disp('Using Iso2Mesh.  Please cite:');
                disp('Qianqian Fang and David Boas, "Tetrahedral mesh generation from volumetric binary')
                disp('and gray-scale images," Proceedings of IEEE International Symposium on Biomedical ');
                disp('Imaging 2009, pp. 1142-1145, 2009');
                
            catch
               
                disp('Please download the iso2mesh package from:');
                disp('http://iso2mesh.sourceforge.net');
                disp('and/or add to the matlab path');
                 error('Cannot find Iso2Mesh on Matlab Path');
                
            end;
            
            [n,f]=v2s(obj.vol,.01*max(obj.vol(:)),1);
            
            %Correct the scaling
            n=n.*(ones(size(n,1),1)*obj.dim);
            n=n+ones(size(n,1),1)*obj.origin;
            
            mesh = nirs.core.Mesh(n,f(:,1:end-1),[]);
            mesh.regions=ones(size(n,1),1);
        end
        
    end
    
end

