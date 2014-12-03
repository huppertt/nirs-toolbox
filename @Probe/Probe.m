classdef Probe
    %PROBE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        srcPos; % nSrc x 3
        srcDir; % nSrc x 3
        detPos; % nDet x 3
        detDir; % nDet x 3
        lambda; % nWavelengths x 1
        link;   % nChannel x [iSrc, iDet, iLambda]
        refPts; % possible set of reference points for alignment w/ image or atlas
        description;
    end
    
    properties( Dependent = true )
        distances;
    end
    
    methods
        %% Constructor
        function obj = Probe( varargin )
            if nargin > 0
                obj.srcPos = varargin{1};
            end
            
            if nargin > 1
                obj.detPos = varargin{2};
            end
            
            if nargin > 2
                obj.link = varargin{3};
            end
            
            if nargin > 3
                obj.lambda = varargin{4};
            end
       
            if nargin > 4
                obj.refPts = varargin{5};
            end
            
            if nargin > 5
                obj.description = varargin{6};
            end
            
            if nargin > 6
                obj.srcDir = varargin{7};
            end
            
            if nargin > 7
                obj.detDir = varargin{8};
            end
            
            if nargin > 8
                error('Too many input arguments.')
            end
        end

        %% Set/Get
        function obj = set.srcPos( obj,newSrcPos )
            if ismatrix(newSrcPos) && size(newSrcPos,2) == 3 || isempty(newSrcPos)
                obj.srcPos = newSrcPos;
            else
                error( 'srcPos should be a n x 3 array.' )
            end
        end
        
        function obj = set.srcDir( obj,newSrcDir )
            if ismatrix(newSrcDir) && size(newSrcDir,2) == 3 || isempty(newSrcDir)
                obj.srcDir = newSrcDir;
            else
                error( 'srcDir should be a n x 3 array.' )
            end
        end
        
        function obj = set.detPos( obj,newDetPos )
            if ismatrix(newDetPos) && size(newDetPos,2) == 3 || isempty(newDetPos)
                obj.detPos = newDetPos;
            else
                error( 'detPos should be a n x 3 array.' )
            end
        end
        
        function obj = set.detDir( obj,newDetDir )
            if ismatrix(newDetDir) && size(newDetDir,2) == 3 || isempty(newDetDir)
                obj.detDir = newDetDir;
            else
                error( 'detDir should be a n x 3 array.' )
            end
        end
        
        function obj = set.link( obj, newLink )
            if ismatrix(newLink) && (size(newLink,2) == 3 || size(newLink,2) == 2 || isempty(newLink))
                obj.link = newLink;
            else
                error('Link should be a matrix with 3 columns [iSrc iDet Lambda].')
            end
        end
        
        function obj = set.lambda( obj, newLambda )
            if isvector( newLambda )
                if iscolumn( newLambda )
                    obj.lambda = newLambda;
                else
                    obj.lambda = newLambda';
                end
            end
        end
        
        function obj = set.description( obj, newDescription )
            if ischar( newDescription ) || isempty(newDescription)
                obj.description = newDescription;
            else
                error( 'Description should be a string.' )
            end
        end 
        
        function obj = set.refPts( obj,newRefPts )
            if ismatrix(newRefPts) && size(newRefPts,2) == 3 || isempty(newRefPts)
                obj.refPts = newRefPts;
            else
                error( 'refPts should be a n x 3 array.' )
            end
        end
        
        %% Methods
        function d = get.distances( obj )
            d = sqrt( sum( (obj.srcPos(obj.link(:,1),:) - obj.detPos(obj.link(:,2),:))'.^2 ))';
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
%             out = out && all( size(obj.srcPos) == size(obj.srcDir) );
        end
        
        function out = isValidDets( obj )
            out = obj.isUniqueDets;
            out = out && (max(obj.link(:,1)) == size(obj.srcPos,1));
            out = out && (max(obj.link(:,2)) == size(obj.detPos,1));
%             out = out && all( size(obj.detPos) == size(obj.detDir) );
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
            
            obj.srcDir = obj.srcDir( uSrc(:,1),: );
            obj.detDir = obj.detDir( uDet(:,1),: );
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

