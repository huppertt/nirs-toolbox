function data = nirstorm2toolbox(sDataIn,events,ChanneMat);

data=nirs.core.Data;
data.data=sDataIn.F';
data.time=linspace(sDataIn.Time(1),sDataIn.Time(2),size(data.data,1));

for i=1:length(events)
    s=nirs.design.StimulusEvents;
    s.name=events(i).label;
    s.onset=events(i).times;
    s.amp=events(i).epochs;
    s.dur=ones(size(s.onset));
    data.stimulus(s.name)=s;
end

for i=1:length( ChanneMat.Channel)
    if strcmp(ChanneMat.Channel(i).Type, 'NIRS')
        ml(i,1)=str2num(ChanneMat.Channel(i).Name(strfind(ChanneMat.Channel(i).Name,'S')+1:...
            strfind(ChanneMat.Channel(i).Name,'D')-1));
        ml(i,2)=str2num(ChanneMat.Channel(i).Name(strfind(ChanneMat.Channel(i).Name,'D')+1:...
            strfind(ChanneMat.Channel(i).Name,'W')-1));
        ml(i,4)=str2num(ChanneMat.Channel(i).Name(strfind(ChanneMat.Channel(i).Name,'WL')+2:end));
        
        SrcPos(ml(i,1),:)=ChanneMat.Channel(i).Loc(:,1)';
        DetPos(ml(i,2),:)=ChanneMat.Channel(i).Loc(:,2)';
    end
    
end
SD.Lambda=ChanneMat.Nirs.Wavelengths';
[~,ml(:,4)]=ismember(ml(:,4),SD.Lambda);

SD.MeasList=ml;
  SD.SrcPos=SrcPos;
  SD.DetPos=DetPos;

data.probe=nirs.util.sd2probe(SD);

  