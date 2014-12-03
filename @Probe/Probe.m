classdef Probe
    %PROBE This object hold nirs probe geometries and reformat source and
    %detector indices for forward models.
    
    properties
        srcPos;         % nSrc x 3
        srcDir;         % nSrc x 3
        
        detPos;         % nDet x 3
        detDir;         % nDet x 3
        
        lambda;         % nWavelengths x 1
        
        link;           % nChannel x [iSrc, iDet, iLambda]
        
        description;
    end
    
    properties( Dependent = true )
        distances;
    end
    
    methods
        %% Constructor
        function obj = Probe( srcPos, detPos, link, lambda )
            if nargin > 0, obj.srcPos = srcPos; end
            if nargin > 1, obj.detPos = detPos; end
            if nargin > 2, obj.link = link; end
            if nargin > 3, obj.lambda = lambda; end
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
        
        % optional vectors descibing the direction of a pencil beam
        function obj = set.srcDir( obj,srcDir )
            assert( ismatrix(srcDir) )
            obj.srcDir = srcDir;
        end
        
        function obj = set.detDir( obj,detDir )
            assert( ismatrix(detDir) )
            obj.detDir = detDir;
        end
        
        % link [iSrc iDet iLambda]
        function obj = set.link( obj,link )
            assert( ismatrix(link) && size(link,2) == 3 )
            obj.link = link;
        end
        
        % wavelengths
        function obj = set.lambda( obj, lambda )
            assert( isvector(lambda) )
            obj.lambda = lambda;
        end
        
        function obj = set.description( obj, description )
            assert( ischar( description ) || isempty( description ) )
          	obj.description = description;
        end
            
        %% Methods
        function d = get.distances( obj )
            isrc = obj.link(:,1);
            idet = obj.link(:,2);
            
            vec = obj.srcPos(isrc,:) - obj.detPos(idet,:);
            
            d = sqrt( sum( vec.^2,2 ) );
        end
        
        function out = isUniqueSrcs( obj )
            % check that sources have unique wavelengths
            [~,~,idx1] = unique(obj.link(:,[1 3]),'rows');
            [~,~,idx2] = unique(obj.link(:,1));
            
            if all( idx1 == idx2 )
                out = 1;
            else
                out = 0;
            end
            
        end
        
        function out = isUniqueDets( obj )
            % check that detectors have unique wavelengths
            [~,~,idx1] = unique(obj.link(:,[2 3]),'rows');
            [~,~,idx2] = unique(obj.link(:,2));
            
            if all( idx1 == idx2 )
                out = 1;
            else
                out = 0;
            end
        end
        
        function out = isUnique( obj )
            out = obj.isUniqueSrcs && obj.isUniqueDets;
        end
        
        function out = isValidSrcs( obj )
            out = obj.isUniqueSrcs;
            out = out && (max(obj.link(:,1)) == size(obj.srcPos,1));
            out = out && (max(obj.link(:,2)) == size(obj.detPos,1));
            
            out = out && (max(obj.link(:,3)) == length(obj.lambda));
        end
        
        function out = isValidDets( obj )
            out = obj.isUniqueDets;
            out = out && (max(obj.link(:,1)) == size(obj.srcPos,1));
            out = out && (max(obj.link(:,2)) == size(obj.detPos,1));
        end
        
        function out = isValid( obj )
            out = obj.isValidSrcs && obj.isValidDets;
        end
        
        function obj = swapSD( obj )
            detPos = obj.srcPos;
            srcPos = obj.detPos;
            detDir = obj.srcDir;
            srcDir = obj.detDir;
            link = obj.link;
            link(:,1) = obj.link(:,2);
            link(:,2) = obj.link(:,1);

            obj.detPos = detPos;
            obj.srcPos = srcPos;
            obj.detDir = detDir;
            obj.srcDir = srcDir;
            obj.link = link; 
        end
        
        function obj = makeUniqueProbe( obj )
            % generates an equivalent probe with unique SD idx's
            [uSrc,~,iSrc] = unique( obj.link(:,[1 3]),'rows' );
            [uDet,~,iDet] = unique( obj.link(:,[2 3]),'rows' );
            
            obj.link(:,1:2) = [iSrc iDet];
            
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
        
%         function display( obj, varargin )
%             if nargin == 1
%                 figure, hold on, 
%                 plot( obj.detPos(:,1), obj.detPos(:,2), 'bo','MarkerSize',10,'LineWidth',2 )
%                 plot( obj.srcPos(:,1), obj.srcPos(:,2), 'rx','MarkerSize',14,'LineWidth',2 )
%                 axis normal
%                 legend( 'Detectors', 'Sources' )
%             end
%         end
end

