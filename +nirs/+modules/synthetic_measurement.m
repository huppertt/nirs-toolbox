classdef synthetic_measurement < nirs.modules.AbstractModule
%% Combines ChannelStats data into a common space
% 
    properties
        commonprobe = 'combine'; 
    end
    methods
        function obj = synthetic_measurement( prevJob )
           obj.name = 'synthetic measurement';
           if nargin > 0
               obj.prevJob = prevJob;
           end
        end
        
        function data = runThis( obj, data )
            if(~isa(obj.commonprobe,'nirs.core.Probe') && ...
                    ~isa(obj.commonprobe,'nirs.core.Probe1020'))
                probe=combineprobes(data);
            else
                probe=obj.commonprobe;
            end
            for i = 1:length(data)
                data(i)=synthetic_meas(data(i),probe);
            end
        end
    end
    
end


function probe = combineprobes(data)

ndet=0;
nsrc=0;
link=[];
optodes=[];
for i=1:length(data)
    thislink=data(i).probe.link;
    thislink.source=thislink.source+nsrc;
    thislink.detector=thislink.detector+ndet;
    
    thisoptodes=data(i).probe.optodes;
    for j=1:size(data(i).probe.detPos,1)
        dI=['0000' num2str(j)];
        dI=dI(end-3:end);
        dI2=['0000' num2str(j+ndet)];
        dI2=dI2(end-3:end);
        
        lst=find(ismember(thisoptodes.Name,['Detector-' dI]));
        thisoptodes.Name{lst}=['Detector-' dI2];
    end
    
    for j=1:size(data(i).probe.srcPos,1)
        dI=['0000' num2str(j)];
        dI=dI(end-3:end);
        dI2=['0000' num2str(j+nsrc)];
        dI2=dI2(end-3:end);
        
        lst=find(ismember(thisoptodes.Name,['Source-' dI]));
        thisoptodes.Name{lst}=['Source-' dI2];
    end
    link=[link; thislink];
    optodes=[optodes; thisoptodes];
    ndet=ndet+size(data(i).probe.detPos,1);
    nsrc=nsrc+size(data(i).probe.srcPos,1);
    
end

probe=data(1).probe;
probe.optodes=optodes;
probe.link=link;

end