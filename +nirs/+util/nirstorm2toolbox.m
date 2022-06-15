function data = nirstorm2toolbox(sDataIn,events,ChanneMat)

data=nirs.core.Data;

nirs_idx = strcmp({ChanneMat.Channel.Type},'NIRS') & (sDataIn.ChannelFlag==1)';
aux_idx = find(~strcmp({ChanneMat.Channel.Type},'NIRS') & (sDataIn.ChannelFlag==1)');

data.data=sDataIn.F(nirs_idx,:)';
data.time=linspace(sDataIn.Time(1),sDataIn.Time(end),size(data.data,1));

for i=1:length(events)
    s=nirs.design.StimulusEvents;
    s.name=events(i).label;
    s.onset=events(i).times(1,:);
    s.amp=events(i).epochs;
    
    if size(events(i).times,1) == 2
        s.dur=events(i).times(2,:)  - events(i).times(1,:);
    else
        s.dur = ones( 1,size(events(i).times,2));
    end
    data.stimulus(s.name)=s;
end

k = 1;
for i=1:length( ChanneMat.Channel)
    if nirs_idx(i)
        ml(k,1)=str2num(ChanneMat.Channel(i).Name(strfind(ChanneMat.Channel(i).Name,'S')+1:...
            strfind(ChanneMat.Channel(i).Name,'D')-1));
        ml(k,2)=str2num(ChanneMat.Channel(i).Name(strfind(ChanneMat.Channel(i).Name,'D')+1:...
            strfind(ChanneMat.Channel(i).Name,'W')-1));
        ml(k,4)=str2num(ChanneMat.Channel(i).Name(strfind(ChanneMat.Channel(i).Name,'WL')+2:end));
        
        SrcPos(ml(k,1),:)=ChanneMat.Channel(i).Loc(:,1)';
        DetPos(ml(k,2),:)=ChanneMat.Channel(i).Loc(:,2)';

        k = k +1 ;
    end
    
end

if ~isempty(aux_idx)
    for j=1:length(aux_idx)
        st=nirs.core.GenericData;
        st.data=sDataIn.F(aux_idx(j),:);
        st.time=data.time;
        
        data.auxillary(ChanneMat.Channel(aux_idx(j)).Name)=st;
    end
end


SD.Lambda=ChanneMat.Nirs.Wavelengths';
[~,ml(:,4)]=ismember(ml(:,4),SD.Lambda);

SD.MeasList=ml;
SD.SrcPos=SrcPos;
SD.DetPos=DetPos;

data.probe=nirs.util.sd2probe(SD);

end