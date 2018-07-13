function raw = loadOxy3(filename)
% This function reads an Artinis *.oxy3 file
% This function is based on the oxysoft2matlab_v1.09 toolbox provided by
% Artinis which is located in the nirs-toolbox/external folder

nirs_data = oxysoft2matlab(filename, 'rawOD',[], false);

SrcPos=nirs_data.transPos';
SrcPos(:,3)=0;
DetPos=nirs_data.receiPos';
DetPos(:,3)=0;

% Create the link table
for i=1:length(nirs_data.ODlabel)
    
    sIdx=nirs_data.ODlabel{i}(strfind(nirs_data.ODlabel{i},'Tx'):strfind(nirs_data.ODlabel{i},'@')-1);
    source(i,1)=find(ismember(upper(nirs_data.TxLabel),upper(sIdx)));
    
    dIdx=nirs_data.ODlabel{i}(strfind(nirs_data.ODlabel{i},'Rx'):strfind(nirs_data.ODlabel{i},'-')-1);
    
    detector(i,1)=find(ismember(upper(nirs_data.RxLabel),upper(dIdx)));
    
    type(i,1)=str2num(nirs_data.ODlabel{i}(strfind(nirs_data.ODlabel{1},'@')+1:strfind(nirs_data.ODlabel{1},'nm')-1));
end
link=table(source,detector,type);  

% Create the probe object
probe=nirs.core.Probe(SrcPos,DetPos,link);

raw=nirs.core.Data;

raw.probe=probe;
raw.data=exp(-nirs_data.OD);  % Convert to a "raw" signal for consistency to the other systems
raw.time=nirs_data.time;


% Now add stimuls events
if(~isempty(nirs_data.events.names))
    stim=Dictionary;
    for i=1:length(nirs_data.events.names)
        s=nirs.design.StimulusEvents;
        s.name=nirs_data.events.names{i};
        s.onset=nirs_data.events.onsets{i}'/raw.Fs;
        s.dur=ones(size(nirs_data.events.onsets{i}'));  % default to 1s duration
        s.amp=ones(size(nirs_data.events.onsets{i}'));  
        stim(s.name)=s;
    end
    raw.stimulus=stim;
end

raw.description=which(filename);

lst=find(any(raw.data==1,1));
raw.data(lst,:)=[];
raw.time(lst)=[];


end
