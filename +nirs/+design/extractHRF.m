function HRF = extractHRF(Stats,basis,duration,Fs)
% This function returns the impulse response or HRF 
% from a nirs.core.ChannelStats variable  
%
% Inputs:
%     Stats -  a nirs.core.ChannelStats or nirs.core.ImageStats variable
%     basis - the basis function used to compute the stats
%     duration (optional) - convolves the IRF (default =1/Fs)
%     Fs - the sample rate to use (default = 4Hz)

if(nargin<4)
    Fs = 4;
end
if(nargin<3)
    duration=1/Fs;
end

if(length(Stats)>1)
    for idx=1:length(Stats)
         HRF(idx) = extractHRF(Stats(idx),basis,duration,Fs);
    end
    return;
end



conditions=Stats.conditions;
for idx=1:length(conditions)
    str=strsplit(conditions{idx},'_');
    str=strcat(str{1:max(1,length(str)-1)});
    conditions{idx}=str;
end
conditions=unique(conditions);    

lenHRF=20;
t=[0:1/Fs:(duration+lenHRF)*length(conditions)];

stimulus=Dictionary();
for idx=1:length(conditions)
    stim=nirs.design.StimulusEvents(conditions{idx},(idx-1)*(duration+lenHRF)+1,duration,1);
    stimulus(conditions{idx})=stim;
end

[X, names] = nirs.design.createDesignMatrix( stimulus, t, basis);
tbl=Stats.table;
tbl=sortrows(tbl,{'source','detector','type','cond'});

beta=reshape(tbl.beta,size(X,2),[]);
Hbeta=X*beta;
tstat=reshape(tbl.tstat,size(X,2),[]);
Htstat=X*tstat;


HRF=nirs.core.Data();
HRF.description=['HRF from basis: ' Stats.description];
HRF.probe=Stats.probe;
HRF.time=t(1:Fs*(duration+lenHRF));

[~,lst]=sortrows(HRF.probe.link,{'source','detector','type'});

data=[];
link=table;
for idx=1:length(conditions)
    lstT=[(idx-1)*(duration+lenHRF)*Fs+1:(idx)*(duration+lenHRF)*Fs];
    data=[data Hbeta(lstT,lst) Htstat(lstT,lst)];
    type=strcat(HRF.probe.link.type(lst),repmat({[':' conditions{idx}]},height(HRF.probe.link),1));
    type2=strcat(HRF.probe.link.type(lst),repmat({[':' conditions{idx} ':tstat' ]},height(HRF.probe.link),1));
    linktmp=[HRF.probe.link; HRF.probe.link];
    linktmp.type=vertcat(type,type2);
    link=[link; linktmp];
end
HRF.probe.link=link;
HRF.data=data;
HRF.demographics=Stats.demographics;


stimulus=Dictionary();
for idx=1:length(conditions)
    stim=nirs.design.StimulusEvents(conditions{idx},1,duration,1);
    stimulus(conditions{idx})=stim;
end

HRF.stimulus=stimulus;



return
