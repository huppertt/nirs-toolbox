function probe = probe_remove_unconnected(probe)
% This function removes any unconnected Src-Det from a probe.  This often
% resulted from the numbering of the cards on a partial Cw6 system

detIdx=1:size(probe.detPos,1);
srcIdx=1:size(probe.srcPos,1);

unUsedSrc = find(~ismember(srcIdx,probe.link.source));
unUsedDet = find(~ismember(detIdx,probe.link.detector));

lstS=find(ismember(probe.optodes.Type,'Source'));
lstD=find(ismember(probe.optodes.Type,'Detector'));

probe.optodes(lstD(unUsedDet),:)=[];
probe.optodes(lstS(unUsedSrc),:)=[];

detIdx(unUsedDet)=[];
srcIdx(unUsedSrc)=[];
[~,probe.link.detector]=ismember(probe.link.detector,detIdx);
[~,probe.link.source]=ismember(probe.link.source,srcIdx);
