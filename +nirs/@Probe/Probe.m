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
        
%         function out = isUniqueSrcs( obj )
%             % check that sources have unique wavelengths
%             [~,~,idx1] = unique(obj.link(:,[1 3]),'rows');
%             [~,~,idx2] = unique(obj.link(:,1));
%             
%             if all( idx1 == idx2 )
%                 out = 1;
%             else
%                 out = 0;
%             end
%             
%         end
%         
%         function out = isUniqueDets( obj )
%             % check that detectors have unique wavelengths
%             [~,~,idx1] = unique(obj.link(:,[2 3]),'rows');
%             [~,~,idx2] = unique(obj.link(:,2));
%             
%             if all( idx1 == idx2 )
%                 out = 1;
%             else
%                 out = 0;
%             end
%         end
%         
%         function out = isUnique( obj )
%             out = obj.isUniqueSrcs && obj.isUniqueDets;
%         end
%         
%         function out = isValidSrcs( obj )
%             out = obj.isUniqueSrcs;
%             out = out && (max(obj.link(:,1)) == size(obj.srcPos,1));
%             out = out && (max(obj.link(:,2)) == size(obj.detPos,1));
%             
%             out = out && (max(obj.link(:,3)) == length(obj.lambda));
%         end
%         
%         function out = isValidDets( obj )
%             out = obj.isUniqueDets;
%             out = out && (max(obj.link(:,1)) == size(obj.srcPos,1));
%             out = out && (max(obj.link(:,2)) == size(obj.detPos,1));
%         end
%         
%         function out = isValid( obj )
%             out = obj.isValidSrcs && obj.isValidDets;
%         end
        
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
        
        function obj = makeUniqueProbe( obj )
            % generates an equivalent probe with unique SD idx's
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
        
        function draw( obj, values, cmap )
%             if nargin == 1
%             figure, hold on, 
%             plot( obj.detPos(:,1), obj.detPos(:,2), 'bo','MarkerSize',10,'LineWidth',2 )
%             plot( obj.srcPos(:,1), obj.srcPos(:,2), 'rx','MarkerSize',14,'LineWidth',2 )
%             axis normal
%             legend( 'Detectors', 'Sources' )
%             end
            
            if nargin < 2
                cmap = [0.8 0.8 0.8];
            elseif nargin < 3
                cmap = evalc('flipud( cbrewer(''div'',''RdBu'',2001) )');
            end

            if nargin == 1
                values = zeros(size(obj.link.source,1),1);
                z = 0;
                cmap = [0.8 0.8 0.8];
            else
                zmax = max( abs(values) );
                z = linspace(-zmax, zmax, size(cmap,1));
                
            end

            link = unique( [obj.link.source obj.link.detector],'rows' );
            
            s = obj.srcPos;
            d = obj.detPos;
            
            gca, hold on
            for iChan = 1:size(link,1)
                
                iSrc = link(iChan,1);
                iDet = link(iChan,2);
                
                x = [s(iSrc,1) d(iDet,1)]';
                y = [s(iSrc,2) d(iDet,2)]';
                
                h(iChan) = line(x,y,'LineWidth',10)
                text(x(1),y(1),['S' num2str(iSrc)], 'FontSize', 20)
                text(x(2),y(2),['D' num2str(iDet)], 'FontSize', 20)
                
            end
            
            
        end
    
    end
end

