function data = loadHamamatsu_NIRO300(filename)
% This function will load a Hamamatsu NIRO-200 data file.  
% The NIRO-200 file is a comma-deliminated text file with very mimumal information pertaining only to
% pre-calculated HbO2/Hb data

tbl=readtable(filename);

time=tbl.elpsec;
d=nan(length(time),40);

link=struct;
optodes=struct;

lst=find(~isnan(time));

for i=1:20
    
    if(ismember(['O2Hb_' num2str(i)],tbl.Properties.VariableNames))
        
        link.source((i-1)*2+1,1)=i;
        link.detector((i-1)*2+1,1)=i;
        link.type{(i-1)*2+1,1}='hbo';
        link.name{(i-1)*2+1,1}=['NIRO200 Channel-' num2str(i)];
        d(:,(i-1)*2+1)=tbl.(['O2Hb_' num2str(i)]);
        
        link.source((i-1)*2+2,1)=i;
        link.detector((i-1)*2+2,1)=i;
        link.type{(i-1)*2+2,1}='hbr';
        d(:,(i-1)*2+2)=tbl.(['HHb_' num2str(i)]);
        link.name{(i-1)*2+2,1}=['NIRO200 Channel-' num2str(i)];
        
        str=['000' num2str(i)];
        optodes.Name{i,1}=['Source-' str(end-3:end)];
        optodes.X(i,1)=i*10;
        optodes.Y(i,1)=15;
        optodes.Z(i,1)=0;
        optodes.Type{i,1}='Source';
        optodes.Units{i,1}='mm';
        optodes.Name{i+10,1}=['Detector-' str(end-3:end)];
        optodes.X(i+10,1)=i*10;
        optodes.Y(i+10,1)=-15;
        optodes.Z(i+10,1)=0;
        optodes.Type{i+10,1}='Detector';
        optodes.Units{i+10,1}='mm';
    end
    
end

link=struct2table(link);
optodes=struct2table(optodes);


lst2=find(tbl.statcode~=0);
d(lst2,:)=NaN;
d(:,height(link)+1:end)=[];

data=nirs.core.Data;
data.data=d(lst,:);
data.time=time(lst);
data.description=filename;

data.probe=nirs.core.Probe;
data.probe.link=link;
data.probe.optodes=optodes;
