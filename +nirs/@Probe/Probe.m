classdef Probe
    %PROBE This object hold nirs probe geometries
    
    properties
        srcPos;     	% nSrc x 3
        srcDir;      	% nSrc x 3
        
        detPos;     	% nDet x 3
        detDir;       	% nDet x 3
        
        % table describing connections with 
        % source, detector, and type columns        
        link;           
        
        %type = 'wavelength'; % wavelength or chromophore
    end
    
    properties( Dependent = true )
        distances;
    end
    
    methods
        %% Constructor
        function obj = Probe( srcPos, detPos, link )
            if nargin > 0, obj.srcPos = srcPos; end
            if nargin > 1, obj.detPos = detPos; end
            if nargin > 2, obj.link = link; end
        end

        %% Set/Get
        % positions
        function obj = set.srcPos( obj,srcPos )
            assert( ismatrix(srcPos) )
            obj.srcPos = srcPos;
        end
        
        function obj = set.detPos( obj,detPos )
            assert( ismatrix(detPos) )
            obj.detPos = detPos;
        end
        
        % optional vectors descibing orienatation of src lasers
        function obj = set.srcDir( obj,srcDir )
            assert( ismatrix(srcDir) )
            obj.srcDir = srcDir;
        end
        
        % optional vectors descibing orienatation of det
        function obj = set.detDir( obj,detDir )
            assert( ismatrix(detDir) )
            obj.detDir = detDir;
        end
        
        % a table with source, detector, wavelength/chromophore
        function obj = set.link( obj,link )
            assert( istable( link ) && isnumeric(link.source) && isnumeric(link.detector) )
            obj.link = link;
        end
        
        %% Methods
        function d = get.distances( obj )
            isrc = obj.link.source;
            idet = obj.link.detector;
            
            vec = obj.srcPos(isrc,:) - obj.detPos(idet,:);
            
            d = sqrt( sum( vec.^2,2 ) );
        end

        % swap sources with detectors
        function obj = swapSD( obj )
            detPos = obj.srcPos;
            srcPos = obj.detPos;
            detDir = obj.srcDir;
            srcDir = obj.detDir;
            link = obj.link;
            
            link.source = obj.link.detector;
            link.detector = obj.link.source;

            obj.detPos = detPos;
            obj.srcPos = srcPos;
            obj.detDir = detDir;
            obj.srcDir = srcDir;
            obj.link = link; 
        end
        
        % generates an equivalent probe with unique SD idx's
        function obj = makeUniqueProbe( obj )
            
            [uSrc,~,iSrc] = unique( ...
                [obj.link.source obj.link.type], ...
                'rows' );
            [uDet,~,iDet] = unique( ...
                [obj.link.detector obj.link.type], ...
                'rows' );
            
            obj.link(:,1:2) = table(iSrc,iDet);
            
            obj.srcPos = obj.srcPos( uSrc(:,1),: );
            obj.detPos = obj.detPos( uDet(:,1),: );
            
            if ~isempty(obj.srcDir)
                obj.srcDir = obj.srcDir( uSrc(:,1),: );
            end
            if ~isempty(obj.detDir)
                obj.detDir = obj.detDir( uDet(:,1),: );
            end
        end
            
    end
end

