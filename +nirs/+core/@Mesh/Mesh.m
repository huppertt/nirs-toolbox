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
        fiducials=table(cell(0,1),cell(0,1),cell(0,1),cell(0,1),cell(0,1),cell(0,1),cell(0,1),...
            'VariableNames',{'Name','X','Y','Z','Type','Units','Draw'});
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
                obj.nodes=node;
                obj.faces=face;
        else
                vol=elemvolume(obj.nodes,obj.elems);
                [node,elem,face]=s2m(obj.nodes,obj.faces,fract,max(vol));
                obj.nodes=node;
                obj.faces=face(:,1:end-1);
                obj.elems=elem(:,1:end-1);
            end
            
        end
        
        function h = draw( obj, values, vmax, thresh, cmap )
            %% draw - displays mesh
            % 
            % Args:
            %     values - (optional) values at each node
            %     vmax   - (optional) max value to display
            %     thresh - (optional) only display values > thresh
            %     cmap   - (optional) colormap
            
            if nargin < 5, cmap     = []; end  
            if nargin < 4, thresh   = []; end
            if nargin < 3, vmax     = []; end
            if nargin < 2, values   = []; end
            
            if(islogical(values));
                values=1*values;
            end
            
          if ~isempty(obj.elems)
               h = nirs.util.plotmesh( obj.nodes, obj.elems, ...
                	values, vmax, thresh, cmap );
            else
                h = nirs.util.plotmesh( obj.nodes, obj.faces, ...
                	values, vmax, thresh, cmap );
          end
            set(h,'FaceAlpha',obj.transparency);
            
            %Draw the fiducial points
            if(height(obj.fiducials)>0)
                hold on;
                lst=find(obj.fiducials.Draw);
                s=scatter3(obj.fiducials.X(lst),obj.fiducials.Y(lst),obj.fiducials.Z(lst),'k','filled');
                hold off;
            end
            
        end
    end
    
end

