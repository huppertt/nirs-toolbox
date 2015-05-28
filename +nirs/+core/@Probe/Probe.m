classdef Probe
    %PROBE This object hold nirs probe geometries
    
    properties
        srcPos	% nSrc x 3
        detPos	% nDet x 3
        
        link	% list of src det connections
	end
    
    properties( Dependent = true )
        distances;
    end
    
    methods
        %% Constructor
        function obj = Probe( srcPos, detPos, link )
            if nargin > 0, obj.srcPos   = srcPos;   end
            if nargin > 1, obj.detPos   = detPos;   end
            if nargin > 2, obj.link     = link;     end
        end

        %% Methods
        % measurement distances in mm
        function d = get.distances( obj )
            isrc = obj.link.source;
            idet = obj.link.detector;
            
            vec = obj.srcPos(isrc,:) - obj.detPos(idet,:);
            
            d = sqrt( sum( vec.^2,2 ) );
        end

        % swap sources with detectors
        function obj = swapSD( obj )
            d = obj.srcPos;
            s = obj.detPos;
            
            obj.srcPos = d;
            obj.detPos = s;
            obj.link(:,[1 2]) = obj.link(:, [2 1]);
        end
    end
end

