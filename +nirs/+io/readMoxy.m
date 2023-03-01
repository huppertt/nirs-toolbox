function data=readMOXY(filename)

folder=fullfile(fileparts(which('nirs.io.readMOXY')),'private');
copyfile(filename,fullfile(folder,'temp.fit'));
system(['java -jar ' folder filesep 'FitCSVTool.jar -b ' folder filesep 'temp.fit ' folder filesep 'fittemp.csv']);
delete(fullfile(folder,'temp.fit'));

fid=fopen(fullfile(folder,'fittemp.csv'),'r');

line=fgetl(fid);
while(~(contains(line,'Data') & contains(line,'record') & contains(line,'timestamp')));
line=fgetl(fid);
end
ncol=length(strfind(line,','));
str='%s';
for i=1:ncol
    str=[str '%s'];
end
c=textscan(line,str,'delimiter',',');
clear C;
for i=1:ncol
    C{1,i}=c{i}(1);
end
cnt=2;
while(1)
    line=fgetl(fid);
    if(isempty(line) | line==-1)
        break;
    end
    c=textscan(line,str,'delimiter',',');
    for i=1:ncol
    C{cnt,i}=c{i}(1);
    end
    cnt=cnt+1;
end
fclose(fid);
delete(fullfile(folder,'fittemp.csv'));
C=cell2table(C);
C=C(ismember(C.C1,{'Data'}) & ismember(C.C3,{'lap','record'}),:);

Cr=C(ismember(C.C3,{'record'}),:);
Cl=C(ismember(C.C3,{'lap'}),:);

Cr=table2cell(Cr);
flds=unique(reshape(Cr(:,4:3:end),[],1));
CC=cell(size(Cr,1),length(flds));
lst2=4:3:size(Cr,2);
for i=1:size(Cr)
    CC(i,:)=repmat({'"-9999"'},1,length(flds));
    [ii,lst]=ismember(flds,{Cr{i,lst2}});
    CC(i,ii)=Cr(i,lst2(lst(ii))+1);
end

for i=1:length(flds)
        if(~isempty(strfind(flds{i},'timestamp')))
            st.timestamp=          cellfun(@(x)str2num(x(2:end-1)),{CC{:,i}}');
        elseif(~isempty(strfind(flds{i},'distance')))
            st.distance=           cellfun(@(x)str2num(x(2:end-1)),{CC{:,i}}');
        elseif(~isempty(strfind(flds{i},'enhanced_speed')))
            st.enhanced_speed=     cellfun(@(x)str2num(x(2:end-1)),{CC{:,i}}');
        elseif(~isempty(strfind(flds{i},'fractional_cadence')))
            st.fractional_cadence= cellfun(@(x)str2num(x(2:end-1)),{CC{:,i}}');
        elseif(~isempty(strfind(flds{i},'cadence')))
            st.cadence=            cellfun(@(x)str2num(x(2:end-1)),{CC{:,i}}');
        elseif(~isempty(strfind(flds{i},'SmO2 Sensor')))
            st.SmO2=               cellfun(@(x)str2num(x(2:end-1)),{CC{:,i}}');
        elseif(~isempty(strfind(flds{i},'THb Sensor')))
            st.THb=               cellfun(@(x)str2num(x(2:end-1)),{CC{:,i}}');
        end
end

if(~isfield(st,'THb'))
    st.THb=nan(size(st.timestamp));
end

if(~isfield(st,'SmO2'))
    st.SmO2=nan(size(st.timestamp));
end

time=st.timestamp-st.timestamp(1);
SmO2=st.SmO2;
THb=st.THb;
SmO2(SmO2==0)=NaN;
SmO2(SmO2==-9999)=NaN;
THb(THb==0)=NaN;
THb(THb==-9999)=NaN;

data=nirs.core.Data;
data.time=time;
data.data=[SmO2 THb];
data.description=filename;
data.probe=nirs.core.ProbeROI({'Leg:SmO2','Leg:THb'});

stim=nirs.design.StimulusEvents;
stim.onset=cellfun(@(x)str2num(x(2:end-1)),Cl.C5)-st.timestamp(1);
stim.onset(find(stim.onset>time(end)))=[];
stim.dur=ones(size(stim.onset));
stim.amp=ones(size(stim.onset));

data.stimulus('lapmarks')=stim;



