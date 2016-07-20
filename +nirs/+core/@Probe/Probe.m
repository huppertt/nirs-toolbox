classdef Probe
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
	end
    
    properties( Dependent = true )
        distances   % (dependent) returns measurement distances
        srcPos      % nsrc x 3 array of source positions (mm)
        detPos      % ndet x 3 array of detector positions (mm)
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
            
            Names=cell(0,1);
            Pos=zeros(0,3);
            Type=cell(0,1);
            Units=cell(0,1);
          
            if nargin > 0; 
                % Sort the sources
                for idx=1:size(srcPos,1)
                   str=['000' num2str(idx)]; 
                   Names{end+1,1}=['Source-' str(end-3:end)];
                    Type{end+1,1}='Source';
                    Pos(end+1,1)=srcPos(idx,1);
                    Pos(end,2)=srcPos(idx,2);
                    Pos(end,3)=srcPos(idx,3);
                    Units{end+1,1}='mm';
                   
                end
            end
            if nargin > 1; 
                % Sort the detectors
                for idx=1:size(detPos,1)
                    str=['000' num2str(idx)]; 
                    Names{end+1,1}=['Detector-' str(end-3:end)];
                    Type{end+1,1}='Detector';
                    Pos(end+1,1)=detPos(idx,1);
                    Pos(end,2)=detPos(idx,2);
                    Pos(end,3)=detPos(idx,3);
                    Units{end+1,1}='mm';
                   
                end
            end
            if nargin > 2, obj.link     = link;     end
            
            obj.optodes=table(Names,Pos(:,1),Pos(:,2),Pos(:,3),Type,Units,...
                'VariableNames',{'Name','X','Y','Z','Type','Units'});
            
        end

        function srcPos = get.srcPos(obj)
            %% This function returns the src pos (in mm)
            tbl=sortrows(obj.optodes,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Source'));
            
            
            srcPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];
            
             %Convert to mm if needed
            lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            srcPos(lstCM,:)=srcPos(lstCM,:)*10;
        end
        
        function detPos = get.detPos(obj)
            %% This function returns the det pos (in mm)
             tbl=sortrows(obj.optodes,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Detector'));
            
            
            
            detPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];
            
            %Convert to mm if needed
            lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            detPos(lstCM,:)=detPos(lstCM,:)*10;
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

            lstD=find(ismember(obj.optodes.Type,'Detector'));
            lstS=find(ismember(obj.optodes.Type,'Source'));
            obj.optodes.Type(lstS)=repmat({'Detector'},length(lstS),1);
            obj.optodes.Type(lstD)=repmat({'Source'},length(lstD),1);
            
            
            obj.link(:,[1 2]) = obj.link(:, [2 1]);
        end
    end
end

