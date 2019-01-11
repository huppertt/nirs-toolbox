function [ContVect,R,names,namesOld,types]=ROIhelper(data,R,names)


if(~iscell(R)); R={R}; end;

if(nargin<3 && ~exist('names'))
    for i=1:length(R)
        names{i}=num2str(i);
    end
end
if ischar( names )
    names = {names};
end

if(isa(data,'nirs.core.Probe') | isa(data,'nirs.core.Probe1020'))
    d=nirs.core.Data;
    d.probe=data;
    d.data=zeros(0,height(data.link));
    data=d;
end


data=data.sorted;

% sort probe
link = data(1).probe.link;


[link, ilink] = sortrows(link, {'type','source', 'detector'});

allSrc=unique(link.source);
allDet=unique(link.detector);

Rall={}; namesAll={};
cnt=1;
for idx=1:length(R)
 
    if(~ismember('weight',R{idx}.Properties.VariableNames))    
        R{idx}.weight=ones(height(R{idx}),1);
        R{idx}.weight=R{idx}.weight./sum(R{idx}.weight);
    end
    if(~ismember('Name',R{idx}.Properties.VariableNames) & ismember('name',R{idx}.Properties.VariableNames))    
        R{idx}.Name=R{idx}.name;
        R{idx}.name=[];
    end
    if(~ismember('Name',R{idx}.Properties.VariableNames))    
        R{idx}.Name=repmat({names{idx}},height(R{idx}),1);
    end
       % Deal with any NaN's    
    if(any(all([isnan(R{idx}.detector) isnan(R{idx}.source)],2)))
        [s,d]=meshgrid(allSrc,allDet);
        n=repmat(cellstr(R{idx}.Name{1}),size(s(:)));
        w=repmat(R{idx}.weight(1),size(s(:)));
        
        R{idx}=[R{idx}; table(s(:),d(:),w(:),n,'VariableNames',{'source','detector','weight','Name'})];
    end
    
    
    
    src=R{idx}.source(isnan(R{idx}.detector));
    weight=R{idx}.weight(isnan(R{idx}.detector));
    Name=R{idx}.Name(isnan(R{idx}.detector));
    
    if(length(src)>0)
        R{idx}=[R{idx}; table(kron(src,ones(length(allDet),1)),...
            kron(ones(length(src),1),allDet),kron(weight,ones(length(allDet),1)),...
            reshape(repmat(Name',length(allDet),1),[],1),...
            'VariableNames',{'source','detector','weight','Name'})];
    end
    
    det=R{idx}.detector(isnan(R{idx}.source));
    weight=R{idx}.weight(isnan(R{idx}.source));
    Name=R{idx}.Name(isnan(R{idx}.source));
    if(length(det)>0)
      
        
        R{idx}=[R{idx}; table(kron(ones(length(det),1),allSrc),...
            kron(ones(length(allSrc),1),det),kron(ones(length(allSrc),1),weight),...
            reshape(repmat(Name',length(allSrc),1),[],1),...
            'VariableNames',{'source','detector','weight','Name'})];
    end
    
    R{idx}((isnan(R{idx}.source) | isnan(R{idx}.detector)),:)=[];
    

    % Table contains the Names for the ROI.
    names2=unique(R{idx}.Name);
    sIdx=find(ismember(R{idx}.Properties.VariableNames,'source'));
    dIdx=find(ismember(R{idx}.Properties.VariableNames,'detector'));
    wIdx=find(ismember(R{idx}.Properties.VariableNames,'weight'));
    
    
    lst=~ismember([R{idx}.source R{idx}.detector],[link.source link.detector],'rows');
    R{idx}(lst,:)=[];
    R{idx}=unique(R{idx});
    
    for idx2=1:length(names2)
        lst=ismember(R{idx}.Name,names2{idx2});
        Rall{cnt}=R{idx}(lst,[sIdx dIdx]);
        namesAll{cnt}=names2{idx2};
        if(~isempty(wIdx))
            ContVect{cnt}=R{idx}.weight(lst);
        end
        cnt=cnt+1;
    end
   
end
R=Rall;
names=namesAll;

if(length(data)>1)
    if(isa(data,'nirs.core.Data'))
        tbl=nirs.core.Data;
        tbl(1)=[];
        for idx=1:length(data)
            tbl=[tbl; nirs.util.roiAverage(data(idx),R,names,splitrois)'];
        end
    else
        tbl=[];
        for idx=1:length(data)
            thistbl=nirs.util.roiAverage(data(idx),R,names,splitrois);
            description = {data(idx).description};
            fileIdx=idx;
            thistbl = [repmat(table(fileIdx,description),height(thistbl),1) thistbl];
            tbl = [tbl; thistbl];
        end
    end
    return
end





    % The region definition is a table, parse it to the contrast vector
    types=unique(link.type);
    
    RNew=cell(length(R)*length(types),1);
    NamesNew=cell(length(R)*length(types),1);
    cnt=1;
     for idx2=1:length(R)
    for idx=1:length(types)
       
            if(iscell(types))
                RNew{cnt}=ismember(link,[R{idx2} table(repmat({types{idx}},height(R{idx2}),1),...
                    'VariableNames',{'type'})]);
                NamesNew{cnt}=[names{idx2} ':' types{idx}];
            else
                RNew{cnt}=ismember(link,[R{idx2} table(repmat([types(idx)],height(R{idx2}),1),...
                    'VariableNames',{'type'})]);
                NamesNew{cnt}=[names{idx2} ':' num2str(types(idx))];
            end
            
            if(exist('ContVect'))
                ContVectNew{cnt}=ContVect{idx2};
            end
            cnt=cnt+1;
        end
    end
    R=RNew;
    namesOld=names;
    names=NamesNew;
    ContVect=ContVectNew;