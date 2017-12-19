function s = roi_stats(obj,region)
% this function returns a ROI stats variable for a spatial ROI average from
% an image

if(iscellstr(region))
    mesh=obj.mesh;
    mask=zeros(length(mesh.nodes),length(region));
    for i=1:length(region)
        for j=1:mesh.labels.count
            if(ismember(region{i},mesh.labels.values{j}.Label))
                
                ii=find(ismember(mesh.labels.values{j}.Label,region{i}));
                
                mask(mesh.labels.values{j}.VertexIndex{ii},i)=1;
                
            end
        end
    end
else
    mask=region;
    clear region;
    for i=1:size(mask,2)
        region{i}=['ROI-' num2str(i)];
    end
    
end

cond=obj.conditions;
types=unique(obj.variables.type);


ball=[];
call=[];

tbl=table;
for i=1:length(cond)
    for j=1:length(types)
        lst=find(ismember(obj.variables.cond,cond{i}) & ismember(obj.variables.type,types{j}));
            c=mask'*obj.covb_chol(lst,:);
            b=mask'*obj.beta(lst);
        ball=[ball; b];
        call=blkdiag(call,c*c');
        tbl=[tbl; table(region,repmat({types{j}},length(region),1),...
            repmat({cond{i}},length(region),1),'VariableNames',{'Region','type','cond'})];
    end
end

s=nirs.core.ChannelStats;
s.description=obj.description;
s.beta=ball;
s.covb=call;
s.variables=tbl;
s.demographics=obj.demographics;
s.dfe=obj.dfe;

s.probe=nirs.core.ProbeROI;
s.probe.RegionNames=region;



