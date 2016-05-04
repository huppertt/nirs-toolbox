function HRF = extractHRF(Stats,basis,duration,Fs,type)
% This function returns the impulse response or HRF 
% from a nirs.core.ChannelStats variable  
%
% Inputs:
%     Stats -  a nirs.core.ChannelStats or nirs.core.ImageStats variable
%     basis - the basis function used to compute the stats
%     duration (optional) - convolves the IRF (default =1/Fs)
%     Fs - the sample rate to use (default = 4Hz)

if(nargin<5)
    type={'hrf'};
end

if(~iscell(type))
    type={type};
end

if(nargin<4)
    Fs = 4;
end
if(nargin<3 | isempty(duration))
    duration=1/Fs;
end

if(~isa(basis,'Dictionary'))
    b=Dictionary();
    b('default')=basis;
    basis=b;
end

if(length(Stats)>1)
    for idx=1:length(Stats)
         HRF(idx) = nirs.design.extractHRF(Stats(idx),basis,duration,Fs,type);
    end
    return;
end



conditions=Stats.conditions;
for idx=1:length(conditions)
    l=strfind(conditions{idx},':');
    if(isempty(l))
        l=length(conditions{idx});
    else
        l=l-1;
    end
    
    
    str=strsplit(conditions{idx}(1:l),'_');
    str=strcat(str{1:max(1,length(str)-1)});
    conditions{idx}=[str conditions{idx}(l+1:end)];
end
conditions=unique(conditions);    

lenHRF=90;
t=[0:1/Fs:(duration+lenHRF)*length(conditions)];

stimulus=Dictionary();
for idx=1:length(conditions)
    stim=nirs.design.StimulusEvents(conditions{idx},(idx-1)*(duration+lenHRF)+1,duration,1);
    stimulus(conditions{idx})=stim;
end

[X, names] = nirs.design.createDesignMatrix( stimulus, t, basis);
tbl=Stats.table;
tbl=sortrows(tbl,{'source','detector','type','cond'});

cnames=tbl.cond(1:size(X,2));
for i=1:length(cnames)
     l=strfind(cnames{i},':');
    if(isempty(l))
        l=length(cnames{i});
    else
        l=l-1;
    end
    
    
    str=strsplit(cnames{i}(1:l),'_');
    num=str{end};
    str=strcat(str{1:max(1,length(str)-1)});
    cnames{i}=[str cnames{i}(l+1:end) '_' num];
end


[i,lst]=ismember(names,cnames);

beta=reshape(tbl.beta,size(X,2),[]);
Hbeta=X*beta(lst,:);
tstat=reshape(tbl.tstat,size(X,2),[]);
Htstat=X*tstat(lst,:);

% Cut off all the zeros at the end
[i,~]=find(X~=0);
i=mod(i,(duration+lenHRF)*Fs);
npts = fix(min(max(i)+10,(duration+lenHRF)*Fs));

HRF=nirs.core.Data();
HRF.description=['HRF from basis: ' Stats.description];
HRF.probe=Stats.probe;
HRF.time=t(1:npts);

[~,lst]=sortrows(HRF.probe.link,{'source','detector','type'});

data=[];
link=table;
for idx=1:length(conditions)
    lstT=fix([(idx-1)*(duration+lenHRF)*Fs+[1:npts]]);
    typ={};
    typ2={};
    if(ismember(lower(type),'hrf'))
        data=[data Hbeta(lstT,lst) ];
         typ=strcat(HRF.probe.link.type(lst),repmat({[':' conditions{idx}]},height(HRF.probe.link),1));
    end
    if(ismember(lower(type),'tstat'))
        data=[data  Htstat(lstT,lst)];
        typ2=strcat(HRF.probe.link.type(lst),repmat({[':' conditions{idx} ':tstat' ]},height(HRF.probe.link),1));
   
    end
        
    linktmp=repmat(HRF.probe.link,length(type),1);
    linktmp.type=vertcat(typ,typ2);
    link=[link; linktmp];
end

n=height(tbl)/length(HRF.probe.link.type)/length(conditions);
link=tbl(1:n:end,1:4);
for i=1:height(link)
    link.cond{i}=link.cond{i}(1:max(strfind(link.cond{i},'_'))-1);
end
link.type=strcat(link.type,repmat({':'},height(link),1),link.cond);
link.cond=[];

HRF.probe.link=link;
HRF.data=data;

if(isempty(Stats.demographics))
    Stats.demographics=Dictionary();
end
HRF.demographics=Stats.demographics;


stimulus=Dictionary();
for idx=1:length(conditions)
    stim=nirs.design.StimulusEvents(conditions{idx},1,duration,1);
    stimulus(conditions{idx})=stim;
end

HRF.stimulus=stimulus;

return
