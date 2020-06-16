classdef ProbeHyperscan
    %% ProbeHyperscan.  This is a class derived from the core Probe class which adds two-person drawing
    
    %% PROBE - This object hold nirs probe geometries.
    %
    % Properties:
    %     optodes-  a table containing the position, name, and type of any
    %     probe points including any reference or anchor points used in
    %     registration
    %     srcPos - (dependent) nsrc x 3 array of source positions (mm)
    %     detPos - (dependent) ndet x 3 array of detector positions (mm)
    %     link   - a table containing the columns 'source', 'detector',
    %              and 'type' describing the the connections between sources
    %              and detectors
    %     distances - (dependent) returns measurement distances
    %
    %  Methods:
    %     swapSD - returns a Probe object with sources exchanged for detectors
    %     draw   - displays a visualization of the probe geometry
    
    properties
        optodes    % table describing the src/det and any additional probe points
        link        % table describing the connections of source/detector pairs
        SubjectLabels;
        RotateMatrix;
    end
    
    properties( Dependent = true )
        distances   % (dependent) returns measurement distances
        srcPos      % nsrc x 3 array of source positions (mm)
        detPos      % ndet x 3 array of detector positions (mm)
    end
    
    properties(Hidden = true)
        fixeddistances
    end
    
    methods
        
        function obj = ProbeHyperscan(probe,SubjectLabels,RotateMatrix)
            
            obj.link=probe.link;
            obj.optodes=probe.optodes;
            obj.fixeddistances=probe.fixeddistances;
            
            if(nargin>1)
                obj.SubjectLabels=SubjectLabels;
            end
            if(nargin>2)
                obj.RotateMatrix=RotateMartix;
            else
                obj.RotateMatrix=cell(length(unique(SubjectLabels)),1);
                y=.5*(max(obj.optodes.Y)-min(obj.optodes.Y));
                
                for i=1:length(obj.RotateMatrix)
                    a=(i-1)*2*pi/length(obj.RotateMatrix);
                    obj.RotateMatrix{i}=[1 0 0 0; 0 1 0 y; 0 0 1 0; 0 0 0 1]'*[cos(a) -sin(a) 0 0; sin(a) cos(a) 0 0; 0 0 1 0; 0 0 0 1];
                end
                
            end
        end
        
        
        function srcPos = get.srcPos(obj)
            %% This function returns the src pos (in mm)
            [tbl,l]=sortrows(obj.optodes,{'Type','Name'});
            
            lst=find(ismember(tbl.Type,'Source'));
            
            
            srcPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];
            srcPos(:,4)=1;
            uSubjs=unique(obj.SubjectLabels);
            for i=1:length(uSubjs)
                lst2=find(ismember(obj.SubjectLabels(l(lst)),uSubjs{i}));
                srcPos(lst2,:)=srcPos(lst2,:)*obj.RotateMatrix{i};
            end
            srcPos(:,4)=[];
            
            
            %Convert to mm if needed
            lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            srcPos(lstCM,:)=srcPos(lstCM,:)*10;
        end
        
        function detPos = get.detPos(obj)
            %% This function returns the det pos (in mm)
            [tbl,l]=sortrows(obj.optodes,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Detector'));
            
            
            
            detPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];
            detPos(:,4)=1;
            uSubjs=unique(obj.SubjectLabels);
            for i=1:length(uSubjs)
                lst2=find(ismember(obj.SubjectLabels(l(lst)),uSubjs{i}));
                detPos(lst2,:)=detPos(lst2,:)*obj.RotateMatrix{i};
            end
            detPos(:,4)=[];
            
            %Convert to mm if needed
            lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            detPos(lstCM,:)=detPos(lstCM,:)*10;
        end
        
        
        
        function varargout=draw( obj, varargin)
            pp=nirs.core.Probe;
            pp.link=obj.link;
            pp.optodes=obj.optodes;
            
            uSubjs=unique(obj.SubjectLabels);
            for i=1:length(uSubjs)
                lst=find(ismember(obj.SubjectLabels,uSubjs{i}));
                xyz=[obj.optodes.X(lst) obj.optodes.Y(lst) obj.optodes.Z(lst)];
                xyz(:,4)=1;
                xyz=xyz*obj.RotateMatrix{i};
                pp.optodes.X(lst)=xyz(:,1);
                pp.optodes.Y(lst)=xyz(:,2);
                pp.optodes.Z(lst)=xyz(:,3);
                
            end
            
            if(nargout==0)
                pp.draw(varargin{:});
            elseif(nargout>0)
                varargout{1}=pp.draw(varargin{:});
            end
            
        end
        function d = get.distances( obj )
            %% distances - Calculates measurement distance for each channel.
            
            if(isa(obj,'nirs.core.ProbeROI'))
                d=obj.fixeddistances;
                return
            end
            
            
            isrc = obj.link.source;
            idet = obj.link.detector;
            
            if iscell(isrc)
                for i = 1:length(isrc)
                    vec = obj.srcPos(isrc{i},:) - obj.detPos(idet{i},:);
                    d(i,1) = mean(sqrt(sum(vec.^2,2)));
                end
            else
                vec = obj.srcPos(isrc,:) - obj.detPos(idet,:);
                d = sqrt( sum( vec.^2,2 ) );
            end
            
            if(~isempty(obj.fixeddistances))
                if(any((d-obj.fixeddistances)~=0))
                    warning('probe distances do not match layout');
                end
                d=obj.fixeddistances;
                return
            end
            
            
        end
        
        function obj = swapSD( obj )
            %% swapSD - Swaps sources for detectors and vice versa.
            
            lstD=find(ismember(obj.optodes.Type,'Detector'));
            lstS=find(ismember(obj.optodes.Type,'Source'));
            obj.optodes.Type(lstS)=repmat({'Detector'},length(lstS),1);
            obj.optodes.Type(lstD)=repmat({'Source'},length(lstD),1);
            
            
            obj.link(:,[1 2]) = obj.link(:, [2 1]);
        end
        
        
    end
    
end