classdef ProbeConnections < nirs.core.Probe
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
        connections
    end
    
   methods (Access = protected)
        function link = get_link_local(obj)
            for i=1:height(obj.connections)
                st=obj.link_probe(obj.connections.start(i),:);
                en=obj.link_probe(obj.connections.end(i),:);
                source{i,1}=['S-' num2str(st.source) ':D-' num2str(st.detector)];
                detector{i,1}=['S-' num2str(en.source) ':D-' num2str(en.detector)];
                type{i,1}=obj.connections.type{i};
            end
            link = table(source,detector,type);
        end
         function link = set_link_local(obj,link)
            if(iscell(link.source))
            s=[];
            for i=1:height(link)
                st=link.source{i};
                en=link.detector{i};
                if(~isempty(strfind(link.type{i},'_')))
                    typ=link.type{i}(1:strfind(link.type{i},'_')-1);
                else
                    typ=link.type{i};
                end
                sIdxA=str2num(st(min(strfind(st,'S-')+2:st(strfind(st,':'))-1)));
                dIdxA=str2num(st(strfind(st,':')+3:end));
                sIdxB=str2num(en(min(strfind(en,'S-')+2:st(strfind(en,':'))-1)));
                dIdxB=str2num(en(strfind(en,':')+3:end));
                s.start(i,1)=find(obj.link_probe.source==sIdxA & obj.link_probe.detector==dIdxA & ...
                    ismember(obj.link_probe.type,typ));
                s.end(i,1)=find(obj.link_probe.source==sIdxB & obj.link_probe.detector==dIdxB & ...
                    ismember(obj.link_probe.type,typ));
                
            end
                s.type=link.type;
                obj.connections=struct2table(s);
                link=obj.link_probe;
            else
                link = link;
            end
        end
    end
    methods
        function obj = ProbeConnections(template)
            %% Probe - Creates a probe object.
            % 
            % Args:
            %     srcPos - (optional) nsrc x 3 array of source positions (mm)
            %     detPos - (optional) ndet x 3 array of detector positions (mm)
            %     link   - (optional) a table containing the columns 'source', 'detector', 
            %              and 'type' describing the the connections between sources 
            %              and detectors
            
            obj.optodes=template.optodes;
            obj.fixeddistances=template.fixeddistances;
            
            obj.link_probe=template.link;

            
        end
    
      
        
    end
end

