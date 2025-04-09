classdef Probe1020Connections < nirs.core.Probe1020
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
                typ=link.type{i}(1:strfind(link.type{i},'_')-1);
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
        function obj = Probe1020Connections(template)
            %% Probe - Creates a probe object.
            
            
            obj.optodes_registered=template.optodes_registered;   % Registered version of optodes
            obj.braindepth=template.braindepth;  % depth of brain for modeling
            obj.opticalproperties=template.opticalproperties;  % optical properties for modeling

            obj.defaultdrawfcn=template.defaultdrawfcn;


            obj.AP_distance=template.AP_distance;  % Distance from Oz to Fpz
            obj.LR_distance=template.LR_distance;  % Distance from LPA to RPA
            obj.IS_distance=template.IS_distance;  % Distance from center to Cz
            obj.labels=template.labels;  % labels of 10-20 points
            obj.pts1020=template.pts1020;  % 10-20 points
            obj.mesh=template.mesh;
            obj.fixeddistances=template.fixeddistances;
            obj.link_probe=template.link;
            obj.optodes=template.optodes;
            
        end
    
        function p = convertProbe1020(obj)
            
            p=nirs.core.Probe1020;
            p.optodes_registered=obj.optodes_registered;   % Registered version of optodes
            p.braindepth=obj.braindepth;  % depth of brain for modeling
            p.opticalproperties=obj.opticalproperties;  % optical properties for modeling

            p.defaultdrawfcn=obj.defaultdrawfcn;


            p.AP_distance=obj.AP_distance;  % Distance from Oz to Fpz
            p.LR_distance=obj.LR_distance;  % Distance from LPA to RPA
            p.IS_distance=obj.IS_distance;  % Distance from center to Cz
            p.labels=obj.labels;  % labels of 10-20 points
            p.pts1020=obj.pts1020;  % 10-20 points
            p.mesh=obj.mesh;
            p.fixeddistances=obj.fixeddistances;
            p.link_probe=obj.link_probe;
            p.optodes=obj.optodes;
        end
        
    end
end

