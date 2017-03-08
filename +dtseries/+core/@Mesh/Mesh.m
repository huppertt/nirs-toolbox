classdef Mesh
    %% MESH - This object hold high density time series channel information and the mesh
    % 
    % Properties:
    %     mesh-  a nirs.core.Mesh data type
    %     link   - a table containing the columns 'vertex', 
    %              and 'type' describing the the source of the data
    %  Methods:
    %     draw   - displays a visualization of the probe geometry
    
    properties
        
        mesh    % table describing the src/det and any additional probe points
        link        % table describing the connections of source/detector pairs
	end
    
    methods
        function obj = Mesh( mesh, link )
            %% Mesh- creates mesh object
            % 
            % Args:
            %     Mesh - (optional) nirs.core.Mesh type
             %     link   - (optional) a table containing the columns 'vertex', 
            %              and 'type' describing the the connections between sources 
            %              and detectors
            
            
          
            if nargin > 0; 
               obj.mesh=mesh;
            end
            if nargin > 1; 
                obj.link=link;
            end
           
            
        end

    end
end

