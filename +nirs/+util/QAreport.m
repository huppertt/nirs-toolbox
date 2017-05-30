function tbl = QAreport(data)
% This function will return a QA report on the data

if(length(data)>1)
    tbl=[];
    for i=1:length(data)
        tbl=[tbl; nirs.util.QAreport(data(i))];
    end
    
    if(isa(data,'nirs.core.ChannelStats'));
        x=horzcat(data.tstat);
        for id=1:size(x,2)
            w(:,id)=zscore(x(:,id));
        end
        p=2*tcdf(-abs(w),size(x,2));
        q=nirs.math.fdr(p);
        
         numBadChannels=sum(1.*(q<0.05),1)';
         
        tbl=[tbl table(numBadChannels)];
    end
    return;
end

t=struct;
if(isa(data,'nirs.core.ChannelStats'));
    if(~isempty(data.description))
        t.description=data.description;
    end
    type=unique(data.variables.type);
    
    for i=1:length(type)
        lst=find(ismember(data.variables.type,type{i}));
        t=setfield(t,['MedianAbsT_' type{i}],mad(data.tstat(lst)));
        t=setfield(t,['MaxAbsT_' type{i}],max(abs(data.tstat(lst))));
        t=setfield(t,['MinAbsT_' type{i}],min(abs(data.tstat(lst))));
        
        t=setfield(t,['MedianAbsBeta_' type{i}],mad(data.beta(lst)));
        t=setfield(t,['MaxAbsBeta_' type{i}],max(abs(data.beta(lst))));
        t=setfield(t,['MinAbsBeta_' type{i}],min(abs(data.beta(lst))));
        
    end
    
    

elseif(isa(data,'nirs.core.Data'))
    if(~isempty(data.description))
        t.description=data.description;
    end
    type=unique(data.probe.link.type);
    
    for i=1:length(type)
        if(~iscellstr(type))
            lst=find(ismember(data.probe.link.type,type(i)));
            tt=num2str(type(i));
        else
            lst=find(ismember(data.probe.link.type,type{i}));
            tt=type{i};
        end
        
        d=data.data(:,lst);
        sni=nirs.math.structnoiseindex(d);
        t=setfield(t,['MedianVar_' tt],median(var(d,[],1)));
        t=setfield(t,['MaxVar_' tt],max(var(d,[],1)));
        t=setfield(t,['MinVar_' tt],min(var(d,[],1)));
        
        t=setfield(t,['MedianSNI_' tt],median(sni));
        t=setfield(t,['MaxSNI_' tt],max(sni));
        t=setfield(t,['MinSNI_' tt],min(sni));
        t=setfield(t,['NumBadChannels_' tt],length(find(sni<2)));
    end
end

tbl=struct2table(t);

return