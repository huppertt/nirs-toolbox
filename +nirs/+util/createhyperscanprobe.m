function probe = createhyperscanprobe(probe)
% This function will create a hyperscanning probe by copying the parameters
% from a single probe

Y=probe.optodes.Y;

probe.optodes.Y=probe.optodes.Y-mean(Y);
probeShift=probe;

probeShift.optodes.Y = probeShift.optodes.Y+2.5*max(abs(Y));

s=find(ismember(probeShift.optodes.Type,'Source'));
for i=1:length(s)
    str=['000' num2str(i+size(probe.srcPos,1))];
    probeShift.optodes.Name{s(i)}=['Source-' str(end-3:end)];
end
s=find(ismember(probeShift.optodes.Type,'Detector'));
for i=1:length(s)
    str=['000' num2str(i+size(probe.detPos,1))];
    probeShift.optodes.Name{s(i)}=['Detector-' str(end-3:end)];
end


link=probe.link;
link.source=link.source+size(probe.srcPos,1);
link.detector=link.detector+size(probe.detPos,1);

probeShift.link=link;

probe.link=[probe.link; probeShift.link];
probe.optodes=[probe.optodes; probeShift.optodes];