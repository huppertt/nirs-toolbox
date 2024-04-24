function roitbl = roi_from_tstat(chanstats,cname,type)

if(nargin<2)
    cname=[];
end

if(nargin<3)
    type=chanstats(1).probe.types(1);
end

if(length(chanstats)>1)
    for i=1:length(roitbl)
        roitbl{i}=nirs.util.roi_from_tstat(chanstats(i),cname,type);
    end
    return
end


if(isempty(cname))
    cname=chanstats.conditions;
end

if(~iscellstr(cname))
    cname=cellstr(cname);
end

lst=find(ismember(chanstats.probe.link.type,type));

ro.source=chanstats.probe.link(lst,:).source;
ro.detector=chanstats.probe.link(lst,:).detector;

roitbl=table;
for i=1:length(cname)
    roit=ro;
    roit.weight=chanstats.ttest(cname(i)).tstat(lst);
    roit.name=repmat(cname(i),length(lst),1);
    roitbl=[roitbl; struct2table(roit)];
end


