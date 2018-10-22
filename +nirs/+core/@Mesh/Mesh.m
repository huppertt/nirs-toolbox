classdef Mesh
    %% MESH - Describes a triangular mesh.
    % 
    % Properties: 
    %     nodes   - m x 3 list of node coordinates
    %     elems   - n x 4 list of node indices defining tetrahedral elements
    %     faces   - p x 3 list of node indices defining triangular faces
    %     regions - m x 1 list of region indicies (e.g. tissue identity)
    %     transparency - alpha value for transperency when drawing
    %     fiducials -  A table of fiducial points for registrations and
    %     drawing
    %
    %  Methods:
    %     draw - draws values on the mesh
    
    properties
        nodes   % m x 3 list of node coordinates
        faces   % n x 4 list of node indices defining tetrahedral elements
        elems   % p x 3 list of node indices defining triangular faces
        regions % m x 1 list of region indicies (e.g. tissue identity)
        
        transparency=1;
        fiducials=table(cell(0,1),zeros(0,1),zeros(0,1),zeros(0,1),cell(0,1),cell(0,1),false(0,1),...
            'VariableNames',{'Name','X','Y','Z','Type','Units','Draw'});
        
        labels=Dictionary;  % label field for storing atlas information or images
        
    end
    
    methods
        function obj = Mesh( n, f, e )
            %% Mesh - Creates a mesh object.
            % 
            % Args:
            %     n - (optional) nodes
            %     f - (optional) faces
            %     e - (optional) elements
    
            if nargin > 0, obj.nodes = n; end
            if nargin > 1, obj.faces = f; end
            if nargin > 2, obj.elems = e; end
        end
        
        function obj = addfiducials(obj,tbl)
            % Code to add fiducial markers into the exisiting table
            
            % obj.fiducials.Draw(:)=false;
            obj.fiducials=[obj.fiducials;...
                [tbl table(repmat(true,height(tbl),1),...
                'VariableNames',{'Draw'})]];
        end
        
        function obj = reducemesh(obj,fract)
            %This function reduces a mesh using the iso2mesh tools
            
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
            
            if(~isempty(obj.faces)  & isempty(obj.elems))
                % Surface mesh
                [node,face]=meshresample(obj.nodes,obj.faces,fract);
                if(~isempty(obj.regions))
                    k=dsearchn(obj.nodes,node);
                    obj.regions=obj.regions(k);
                end
                obj.nodes=node;
                obj.faces=face;
            else
                vol=elemvolume(obj.nodes,obj.elems);
                [node,elem,face]=s2m(obj.nodes,obj.faces,fract,max(vol));
                if(~isempty(obj.regions))
                    k=dsearchn(obj.nodes,node);
                    obj.regions=obj.regions(k);
                end
                obj.nodes=node;
                obj.faces=face(:,1:end-1);
                obj.elems=elem(:,1:end-1);
            end
            
        end
        
        function mesh = convert2FEM(obj)
                         % converts a BEM mesh to a FEM mesh
             if~(~isempty(obj(1).faces)  & isempty(obj(1).elems))
                warning('mesh is already in BEM form');
                mesh=obj;
                return;
             end
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
            
            mesh = nirs.core.Mesh;
                
                mesh=obj(1);
                
                
                for idx=2:length(obj)
                    n=size(mesh.nodes,1);
                    mesh.nodes=[mesh.nodes; obj(idx).nodes];
                    mesh.faces=[mesh.faces; obj(idx).faces+n];
                    mesh.elems=[mesh.elems; obj(idx).elems+n];
                    try; mesh.fiducials=[mesh.fiducials; obj(idx).fiducials]; end;
                end
                fid=mesh.fiducials;
                [node,elem,face]=cgals2m(mesh.nodes,obj(1).faces,4,10);
                mesh=nirs.core.Mesh(node(:,1:3),face(:,1:3),elem(:,1:4));
                                
               
                image=convert2image(obj);
                n=node(:,1:3);
                n=round((n.*(ones(size(n,1),1)*image.dim)+ones(size(n,1),1)*image.origin));
                lst=sub2ind(image.size,n(:,1),n(:,2),n(:,3));
                mesh.regions=max(image.vol(lst),1);
                mesh.fiducials=fid;
                
        end
            
        function image = convert2image(obj)
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
            
            
            p0=inf;
            p1=-inf;
            for i=1:length(obj)
                p0=floor(min(min(obj(i).nodes),p0));
                p1=ceil(max(max(obj(i).nodes),p1));
            end
            dx=1;  % 1mm resolution
            
            aseg=surf2vol(obj(1).nodes,obj(1).faces,p0(1)-dx:dx:p1(1)+dx,p0(2)-dx:dx:p1(2)+dx,p0(3)-dx:dx:p1(3)+dx);
            aseg=imfill(aseg,'holes');
            for id=2:length(obj)
                img=surf2vol(obj(id).nodes,obj(id).faces,p0(1)-dx:dx:p1(1)+dx,p0(2)-dx:dx:p1(2)+dx,p0(3)-dx:dx:p1(3)+dx);
                img=id*imfill(img,'holes');
                aseg=max(aseg,img);
            end
            
           image=nirs.core.Image;
           image.vol=aseg;
           image.dim=[dx dx dx];
           image.origin=[find(p0(1)-dx:dx:p1(1)+dx==0) find(p0(2)-dx:dx:p1(2)+dx==0) find(p0(3)-dx:dx:p1(3)+dx==0)];
        end
                
        
        
        function h = draw( obj, values, vmax, thresh, cmap,axis_handle)
            %% draw - displays mesh
            % 
            % Args:
            %     values - (optional) values at each node
            %     vmax   - (optional) max value to display
            %     thresh - (optional) only display values > thresh
            %     cmap   - (optional) colormap
            
            if nargin < 6, axis_handle=gca; end;
            if nargin < 5, cmap     = []; end  
            if nargin < 4, thresh   = []; end
            if nargin < 3, vmax     = []; end
            if nargin < 2, values   = []; end
            
            if(length(obj)>1)
                cnt=0;
                for i=1:length(obj)
                    if(~isempty(values))
                        valuesi=values(cnt+[1:size(obj(i).nodes,1)]);
                    else
                        valuesi=[];
                    end
                    cnt=cnt+size(obj(i).nodes,1);
                    h(i) = draw( obj(i), valuesi, vmax, thresh, cmap,axis_handle );
                    hold(axis_handle,'on');
                end
                return
            end
            
            if(islogical(values));
                values=1*values;
            end
            
          if ~isempty(obj.elems)
               h = nirs.util.plotmesh( obj.nodes, obj.elems, ...
                	values, vmax, thresh, cmap,axis_handle );
            else
                h = nirs.util.plotmesh( obj.nodes, obj.faces, ...
                	values, vmax, thresh, cmap,axis_handle );
          end
            set(h,'FaceAlpha',obj.transparency);
            
            %Draw the fiducial points
            if(height(obj.fiducials)>0)
                hold(axis_handle,'on');
                lst=find(obj.fiducials.Draw);
                s=scatter3(axis_handle,obj.fiducials.X(lst),obj.fiducials.Y(lst),obj.fiducials.Z(lst),'k','filled');
                hold(axis_handle,'off');
            end
            %set(gcf,'color','k');
            set(axis_handle,'color','k');
            
            %rotate3d(axis_handle,'on');
        end
    end
    
end

