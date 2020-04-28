function raw = loadShimadzu(filename,probe)
% Put probe.mat file in root directory

if (nargin<2)
    load probe
end

datatype='raw';
fid = fopen(filename,'r');
line='test';
while(~isnumeric(line) & isempty(strfind(line,'Time(sec)')))
    line=fgetl(fid);
    if(~isempty(strfind(line,'Data Type	Hb')))
        datatype='hb';
    end
end




numTabs=length(find(double(line)==9))+1;
s=''; t='';
for i=1:numTabs
    s=[s '%f'];
    t=[t '%s'];
end
fgetl(fid);
c=textscan(fid,s);

names=textscan(line,t);
for i=5:length(names)
    n{i-4}=names{i}{1};
end

if(strcmp(datatype,'hb'))
    n(ismember(n,'totalHb'))=repmat({'hbt'},length(find(ismember(n,'totalHb'))),1);
    n(ismember(n,'oxyHb'))=repmat({'hbo'},length(find(ismember(n,'oxyHb'))),1);
    n(ismember(n,'deoxyHb'))=repmat({'hb'},length(find(ismember(n,'deoxyHb'))),1);
else
    for i=1:length(n)
        n2(i)=str2num(n{i}(4:strfind(n{i},'nm')-1));
    end
    n=n2;
end

Task=c{2};
Mark=c{3};
Count=c{4};

raw=nirs.core.Data;
raw.description=filename;
raw.probe=probe;
raw.data=horzcat(c{5:end});
raw.time=c{1};

nt=length(unique(n));
% raw.probe.link=table(reshape(repmat(raw.probe.link.source,1,nt),[],1),...
%     reshape(repmat(raw.probe.link.detector,1,nt),[],1),n','VariableNames',{'source','detector','type'});

events=unique(Task);
events(events==0)=[];
for i=1:length(events)
     raw.stimulus(['Task' num2str(events(i))]) = nirs.design.vector2event(raw.time,1*(Task==events(i)),['Task' num2str(events(i))]);
end


events=unique(Mark);
events(events==0)=[];
for i=1:length(events)
     raw.stimulus(['Mark' num2str(events(i))]) = nirs.design.vector2event(raw.time,1*(Mark==events(i)),['Mark' num2str(events(i))]);
end


events=unique(Count);
events(events==0)=[];
for i=1:length(events)
     raw.stimulus(['Count' num2str(events(i))]) = nirs.design.vector2event(raw.time,1*(Count==events(i)),['Count' num2str(events(i))]);
end

raw=raw.sorted();
fclose(fid);

patfile=dir([filename(1:strfind(filename,'.')-1) '*.pat']);
if(~isempty(patfile))
    % read the demogrpahics
    fid=fopen(patfile(1).name,'r');
    line=fgetl(fid);
    while(isempty(strfind(line,'[Section Study]')) & ~isnumeric(line))
        line=fgetl(fid);
        field=line(1:strfind(line,'=')-1);
        
        val=line(strfind(line,'=')+1:end);
        if(~isempty(field))
        raw.demographics(field)=val;
        end
    end
end
    


