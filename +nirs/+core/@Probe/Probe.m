classdef Probe
    %% PROBE - This object hold nirs probe geometries.
    % 
    % Properties: 
    %     srcPos - nsrc x 3 array of source positions (mm)
    %     detPos - ndet x 3 array of detector positions (mm)
    %     link   - a table containing the columns 'source', 'detector', 
    %              and 'type' describing the the connections between sources 
    %              and detectors
    %     distances - (dependent) returns measurement distances
    %     
    %  Methods:
    %     swapSD - returns a Probe object with sources exchanged for detectors
    %     draw   - displays a visualization of the probe geometry
    
    properties
        srcPos      % nsrc x 3 array of source positions (mm)
        detPos      % ndet x 3 array of detector positions (mm)
        
        link        % table describing the connections of source/detector pairs
	end
    
    properties( Dependent = true )
        distances   % (dependent) returns measurement distances
    end
    
    methods
        function obj = Probe( srcPos, detPos, link )
            %% Probe - Creates a probe object.
            % 
            % Args:
            %     srcPos - (optional) nsrc x 3 array of source positions (mm)
            %     detPos - (optional) ndet x 3 array of detector positions (mm)
            %     link   - (optional) a table containing the columns 'source', 'detector', 
            %              and 'type' describing the the connections between sources 
            %              and detectors
                     
            if nargin > 0, obj.srcPos   = srcPos;   end
            if nargin > 1, obj.detPos   = detPos;   end
            if nargin > 2, obj.link     = link;     end
        end

        function d = get.distances( obj )
            %% distances - Calculates measurement distance for each channel.
        
            isrc = obj.link.source;
            idet = obj.link.detector;
            
            vec = obj.srcPos(isrc,:) - obj.detPos(idet,:);
            
            d = sqrt( sum( vec.^2,2 ) );
        end
        
        function obj = swapSD( obj )
            %% swapSD - Swaps sources for detectors and vice versa.

            d = obj.srcPos;
            s = obj.detPos;
            
            obj.srcPos = d;
            obj.detPos = s;
            obj.link(:,[1 2]) = obj.link(:, [2 1]);
        end
    end
end

