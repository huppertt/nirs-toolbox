function C = makeConnectivityTimeCourse(data)
% This function computes the time-dependent correlation between all channels and returns
% it as a core.Data variable, which can then be used in the GLM model
% 
% TODO:
% Support the drawing code and probe methods

if(length(data)>1)
    for i=1:length(data)
        disp([num2str(i) ' of ' num2str(length(data))]);
        C(i)=nirs.util.makeConnectivityTimeCourse(data(i));
    end
    return
end


C=data;

probe=data.probe;
types=probe.types;

link=struct;
optodes=struct;
cnt=1; 

optodes.Name={}; optodes.X=[]; optodes.Y=[]; optodes.Z=[]; 
optodes.Type={}; optodes.Units={}; 
optodes.OriginalSource={}; optodes.OriginalDetector={};


for itype=1:length(types)
    lst=find(ismember(probe.link.type,types(itype)));
    
    for i=1:length(lst)
        for j=i+1:length(lst)
            a=data.data(:,lst(i));
            b=data.data(:,lst(j));
            a=a-mean(a); a=a-mean(a);  % do twice to fix numerical rounding issues
            b=b-mean(b); b=b-mean(b);
            a=a/sqrt(var(a));
            b=b/sqrt(var(b));
            
            d(:,cnt)=a.*b;
            link.source(cnt,1)=i;
            link.detector(cnt,1)=j-1;
            link.type(cnt,1)=types(itype);
            
            str=['000' num2str(i)];
            optodes.Name{end+1,1}=['Source-' str(end-3:end)];
            str=['000' num2str(probe.link.source(lst(i)))];
            optodes.OriginalSource{end+1,1}=['Source-' str(end-3:end)];
            str=['000' num2str(probe.link.detector(lst(i)))];
            optodes.OriginalDetector{end+1,1}=['Detector-' str(end-3:end)];
            
            xyz=.5*(probe.srcPos(probe.link.source(lst(i),:),:)+...
                probe.detPos(probe.link.detector(lst(i),:),:));
            optodes.X(end+1,1)=xyz(1);
            optodes.Y(end+1,1)=xyz(2);
            optodes.Z(end+1,1)=xyz(3);
            optodes.Type{end+1,1}='Source';
            optodes.Units{end+1,1}='mm';
            
            
            str=['000' num2str(j-1)];
            optodes.Name{end+1,1}=['Detector-' str(end-3:end)];
            str=['000' num2str(probe.link.source(lst(j)))];
            optodes.OriginalSource{end+1,1}=['Source-' str(end-3:end)];
            str=['000' num2str(probe.link.detector(lst(j)))];
            optodes.OriginalDetector{end+1,1}=['Detector-' str(end-3:end)];
            
            xyz=.5*(probe.srcPos(probe.link.source(lst(j),:),:)+...
                probe.detPos(probe.link.detector(lst(j),:),:));
            optodes.X(end+1,1)=xyz(1);
            optodes.Y(end+1,1)=xyz(2);
            optodes.Z(end+1,1)=xyz(3);
            optodes.Type{end+1,1}='Detector';
            optodes.Units{end+1,1}='mm';
            
            cnt=cnt+1;
        end
    end

end

C.probe.optodes=unique(struct2table(optodes));
C.probe.link=struct2table(link);
C.data=d;
C=C.sorted;




