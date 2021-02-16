function probe = probe_remove_unconnected(probe)
% This function removes any unconnected Src-Det from a probe.  This often
% resulted from the numbering of the cards on a partial Cw6 system


Src=probe.optodes(ismember(probe.optodes.Type,'Source'),:).Name;
Det=probe.optodes(ismember(probe.optodes.Type,'Detector'),:).Name;

srcIdx=zeros(length(Src),1);
detIdx=zeros(length(Det),1);
for i=1:length(Src)
    srcIdx(i)=str2num(Src{i}(strfind(Src{i},'-')+1:end));
end
for i=1:length(Det)
    detIdx(i)=str2num(Det{i}(strfind(Det{i},'-')+1:end));
end

unUsedSrc = find(~ismember(srcIdx,probe.link.source));
unUsedDet = find(~ismember(detIdx,probe.link.detector));

unUsed={Src{unUsedSrc} Det{unUsedDet}};

probe.optodes(ismember(probe.optodes.Name,unUsed),:)=[];
 
if(isa(probe,'nirs.core.Probe1020'))
  probe.optodes_registered(ismember(probe.optodes_registered.Name,unUsed),:)=[];
end

