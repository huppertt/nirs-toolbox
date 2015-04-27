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
        
        % draw probe with optional channel values
        function draw( obj, values, vmax, thresh, cmap )

            % specify colormap
            if nargin == 1
                cmap = [0.5 0.5 1];
            elseif nargin < 5 || isempty(cmap)
                [~,cmap] = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            end
            
            % no thresholding if not specified
            if nargin < 4
                thresh = 0;
            end

            % calculate the mapping from values to colormap
            if nargin == 1
                values = zeros(size(obj.link.source,1),1);
                z = 0; zmax = 0;
            elseif nargin < 4
                zmax = ceil( max( abs(values) ) );
                z = linspace(-zmax, zmax, size(cmap,1));
            else
             	zmax = vmax;
                z = linspace(-zmax, zmax, size(cmap,1));
            end
            
            % threshold colormap
            lst = abs(z) < thresh;
            [~,i] = min(abs(z));
            cmap(lst,:) = repmat( cmap(i,:), [sum(lst) 1] );

            % loop through channels and draw lines
            link = unique( [obj.link.source obj.link.detector],'rows' );
            
            s = obj.srcPos;
            d = obj.detPos;
            
            gca, hold on
            for iChan = 1:size(link,1)
                
                iSrc = link(iChan,1);
                iDet = link(iChan,2);
                
                x = [s(iSrc,1) d(iDet,1)]';
                y = [s(iSrc,2) d(iDet,2)]';
                
                
                [~, iCol] = min(abs(values(iChan)-z));
                h = line(x,y,'LineWidth',5,'Color',cmap(iCol,:));
                text(x(1),y(1),['S' num2str(iSrc)], 'FontSize', 14)
                text(x(2),y(2),['D' num2str(iDet)], 'FontSize', 14)
            end
            
            axis tight
           	axis off
            
            % colormap/bar
            colormap(cmap);
            if zmax ~= 0
                c = colorbar; 
                caxis([-zmax zmax]);
            end
            
            % adjust axes
            pi = get(gca,'Position');
            po = get(gca,'OuterPosition');
            
            po(3:4) = 1.2*pi(3:4);
            set(gca,'OuterPosition',po);
            
            xl = xlim;
            yl = ylim;

            xl = 1.2*diff(xl)/2*[-1 1]+mean(xl);
            yl = 1.2*diff(yl)/2*[-1 1]+mean(yl);

            axis([xl yl])
                
        end
    
    end
end

