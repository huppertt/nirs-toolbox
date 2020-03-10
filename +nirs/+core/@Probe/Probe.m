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
        types
    end
    
    properties(Hidden = true)
        fixeddistances
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
            
            fixeddistances=[];
            regionnames={};
            isregionofinterest=false;
            
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

        function types = get.types(obj)
            types=unique(obj.link.type);
        end
        
        function srcPos = get.srcPos(obj)
            %% This function returns the src pos (in mm)
            
            optodes=obj.optodes;
            lst=[];
            for i=1:height(optodes)
                if(isempty(optodes.Name{i}))
                    lst=[lst; i];
                end
            end
            optodes(lst,:)=[];
            
            tbl=sortrows(optodes,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Source'));
            
            
            srcPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];
            
             %Convert to mm if needed
            lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            srcPos(lstCM,:)=srcPos(lstCM,:)*10;
        end
        
        function detPos = get.detPos(obj)
            
             optodes=obj.optodes;
            lst=[];
            for i=1:height(optodes)
                if(isempty(optodes.Name{i}))
                    lst=[lst; i];
                end
            end
            optodes(lst,:)=[];
            
            %% This function returns the det pos (in mm)
             tbl=sortrows(optodes,{'Type','Name'});
            lst=find(ismember(tbl.Type,'Detector'));
            
            
            
            detPos=[tbl.X(lst) tbl.Y(lst) tbl.Z(lst)];
            
            %Convert to mm if needed
            lstCM=find(ismember(tbl.Units(lst),{'cm'}));
            detPos(lstCM,:)=detPos(lstCM,:)*10;
        end
        
        function d = get.distances( obj )
            %% distances - Calculates measurement distance for each channel.
        
           if(~isempty(obj.fixeddistances))
               d=obj.fixeddistances;
               if(length(d)<height(obj.link))
                    dd=zeros(height(obj.link),1);
                    t=unique(obj.link.type);
                    for i=1:length(t)
                        dd(ismember(obj.link.type,t(i)))=d;
                    end
                    d=dd;
               end
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
            
%             if(~isempty(obj.fixeddistances))
% %                 if(any((d-obj.fixeddistances)~=0))
% %                     warning('probe distances do not match layout');
% %                 end
%                 d=obj.fixeddistances;
%                 return
%             end
%             
            
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

