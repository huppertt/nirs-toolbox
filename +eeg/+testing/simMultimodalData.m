function [data, truth, truthImage] = simMultimodalData(noise,stim,fwdmodel,beta)
% This function simulates multimodal data (currently NIRS & EEG)

if(nargin<1 || isempty(noise))
    probe=nirs.testing.simProbe;
    data.eeg=eeg.testing.simARNoise;
    data.nirs=nirs.testing.simARNoise(probe);
else
    data.eeg=eeg.core.Data; data.eeg(:)=[];
    data.nirs=nirs.core.Data; data.nirs(:)=[];
    for i=1:length(noise)
        if(isa(noise{i},'nirs.core.Data'))
            data.nirs=noise{i};
        else(isa(noise{i},'eeg.core.Data'))
             data.eeg=noise{i};
        end
    end
    
end

time=[];
if(~isempty(data.nirs)); time=[time; data.nirs.time]; end;
if(~isempty(data.eeg)); time=[time; data.eeg.time]; end;
time=sort(unique(time));

if nargin < 2 || isempty(stim)
    stim = nirs.testing.randStimDesign(time, .1, 7, 1);
end

if(nargin<3 || isempty(fwdmodel))
    mesh=nirs.registration.Colin27.BEM;
    fwd.nirs=nirs.forward.ApproxSlab;
    fwd.nirs.Fm=0;
    fwd.nirs.mesh=mesh.mesh(3);
    lambda=unique(data.nirs.probe.link.type);
    fwd.nirs.prop=nirs.media.tissues.brain(.7,50,lambda);
    fwd.nirs.probe=data.nirs.probe;
    
    fwd.eeg=eeg.forward.FieldTrip;
    fwd.eeg.mesh=mesh.mesh;
    fwd.eeg.mesh(4)=[];
    fwd.eeg.probe=data.eeg.probe;
    fwd.eeg.prop=[1 NaN NaN];

end


if nargin < 4 || isempty(beta)
    if(isfield(fwd,'nirs'))
        mesh=fwd.nirs.mesh;
        
        pos=(fwd.nirs.probe.srcPos(fwd.nirs.probe.link.source,:)+...
            fwd.nirs.probe.detPos(fwd.nirs.probe.link.detector,:))/2;
        
    end
    mesh(end).transparency=1;
    nVox=size(mesh(end).nodes,1);
    beta=zeros(nVox,length(stim.keys));
    
    for idx=1:length(stim.keys)
        p=pos(randi(size(pos,1),1,1),:);
        d=sqrt(sum((mesh(end).nodes-ones(nVox,1)*p).^2,2));
        lst=find(d<20 & d>5);
        beta(lst,idx)=1;
    end
end

[rawNIRS,truthImage,~,truthchanNIRS]=nirs.testing.simDataImage(fwd.nirs,data.nirs,stim,{'BA-46_R'});
truthImage=truthImage(1:end/2);
[rawEEG,~,~,truthchanEEG]=eeg.testing.simDataImage(fwd.eeg,data.eeg,stim,{'BA-46_R'});

data={};
data{1}=rawEEG;
data{2}=rawNIRS;

truth{1}=truthchanEEG;
truth{2}=truthchanNIRS;


