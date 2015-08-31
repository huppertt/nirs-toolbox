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
            
          if isempty(obj.faces)
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

