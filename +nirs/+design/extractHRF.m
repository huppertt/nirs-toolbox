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
    if(~isempty(l))
       conditions{idx}=conditions{idx}(1:l-1);
    end
    
end
conditions=unique(conditions);    

if(isa(duration,'Dictionary'))
    dur=zeros(length(conditions),1);
    for idx=1:length(conditions)
        k=find(ismember(duration.keys,conditions{idx}));
        if(~isempty(k))
        s=duration(duration.keys{k});
        dur(idx)=s.dur;
        else
            dur(idx)=1/Fs;
        end
    end
    duration=dur;
end
if(length(duration)==1)
    duration=repmat(duration,length(conditions),1);
end


lenHRF=90;
t=[0:1/Fs:(max(duration)+lenHRF)*length(conditions)];

stimulus=Dictionary();
for idx=1:length(conditions)
    stim=nirs.design.StimulusEvents(conditions{idx},0,duration(idx),1);
    stimulus(conditions{idx})=stim;
end

[X, names] = nirs.design.createDesignMatrix( stimulus, t, basis);
tbl=Stats.table;
tbl=sortrows(tbl,{'source','detector','type','cond'});

cnames=unique(tbl.cond);
cnamesOrig=cnames;
for i=1:length(cnames)
   l=strfind(cnames{i},':');
    if(~isempty(l))
        l=length(cnames{i});
        
        str=strsplit(cnames{i}(1:l),'_');
        num=str{end};
        if(~isempty(str2num(num)))
            str=strcat(str{1:max(1,length(str)-1)});
            cnames{i}=[str cnames{i}(l+1:end) '_' num];
            
        end
    end
end


for j=1:length(cnames)
    cnames{j}(end+1)='&';   % mark the end so stim1 and stim10 are not mistaken
end
for j=1:length(names)
    names{j}(end+1)='&';
end

cnt=1; n={};
for i=1:length(cnames)
    for j=1:length(names)
        if(~isempty(strfind(cnames{i},names{j})))
            n{cnt}=cnames{i}(length(names{j})+1:end);
            cnt=cnt+1;
           
        end
    end
end
no=n;
n=unique(n);

n2={}; cnt=1;
c={};
X = repmat(X,1,length(n));
for i=1:length(n)
    for j=1:length(names)
        n2{cnt,1}=[names{j}(1:end-1) n{i}];
        cnt=cnt+1;
    end
    for j=1:length(conditions)
        c{end+1}=[conditions{j} n{i}];
    end
end
names=n2;
conditions2=c;
lstC={};
for i=1:length(conditions2)
    
    if(length(n)==1)
        nI=i;
        cI=i;
    else
        nI=mod(i-1,length(conditions))+1;
        cI=floor((i)/(length(n)+1))+1;
    end
     lstC{i}=[];
    for j=1:length(names)
        nn=['000' num2str(j)];
        nn=nn(end-1:end);
        nn=[':' nn];
        for k=1:length(names)
            if(~isempty(no{nI}))
                if(~isempty(strfind(names{k},nn)) & ~isempty(strfind(names{k},no{nI})) &...
                        ~isempty(strfind(names{k},[conditions{cI}])))
                    lstC{i}=[lstC{i} k];
                end
            else
                 if(~isempty(strfind(names{k},nn)) &...
                        ~isempty(strfind(names{k},[conditions{cI} ':'])))
                    lstC{i}=[lstC{i} k];
                 elseif(strcmp(names{k},[conditions{cI}]))
                     lstC{i}=k;
                 end
            end
        end
    end
    lstC{i}=unique(lstC{i});
end

for i=1:length(cnames)
    cnames{i}(end)=[];   % mark the end so stim1 and stim10 are not mistaken
end
   

[i,lst2]=ismember(names,cnames);
l=unique(tbl(:,1:3));

% Cut off all the zeros at the end
[i,~]=find(X~=0);
i=mod(i,(max(duration)+lenHRF)*Fs);
npts = fix(min(max(i)+10,(max(duration)+lenHRF)*Fs));

HRF=nirs.core.Data();
HRF.description=['HRF from basis: ' Stats.description];
HRF.probe=Stats.probe;
HRF.time=t(1:npts);
lstT=[1:npts];

link=table;
data=[];



for i=1:height(l)
    if(isa(l.type,'cell'))
        ltype=l.type{i};
        ltype2=ltype;
    else
        ltype=l.type(i);
        ltype2=num2str(ltype);
    end
     lst=find(tbl.source==l.source(i) & tbl.detector==l.detector(i) & ismember(tbl.type,ltype));
     [~,or]=ismember({cnamesOrig{lst2}},tbl(lst,:).cond);
     for j=1:length(conditions2)
         Hbeta=X(lstT,lstC{j})*tbl.beta(lst(or(lstC{j})),:);
         Htstat=X(lstT,lstC{j})*tbl.tstat(lst(or(lstC{j})),:);
         if(ismember(lower(type),'hrf'))
             data=[data Hbeta ];
             
         end
         if(ismember(lower(type),'tstat'))
             data=[data  Htstat];
             
         end
         link=[link; table(l.source(i),l.detector(i),cellstr([ltype2 '_' conditions2{j}]),'VariableNames',{'source','detector','type'})]; 
         
     end
end

HRF.probe.link=link;
HRF.data=data;

if(isempty(Stats.demographics))
    Stats.demographics=Dictionary();
end
HRF.demographics=Stats.demographics;


stimulus=Dictionary();
for idx=1:length(conditions)
    stim=nirs.design.StimulusEvents(conditions{idx},.001,duration(idx),1);
    stimulus(conditions{idx})=stim;
end

HRF.stimulus=stimulus;

return