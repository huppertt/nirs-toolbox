function Stats = add_type_interactions(Stats)
% This function changes the link.type and condition names to allow Ttests
% between different types (e.g. Hbo vs HbR).  
%   
% Stats = nirs.util.add_type_interactions(Stats);
% Stats.ttest('A_hbo+A_hbr');


if(length(Stats)>1)
    for i=1:length(Stats)
         Stats(i) = nirs.util.add_type_interactions(Stats(i));
    end
    return
end

types=Stats.probe.types;

vars=Stats.variables;
for i=1:height(vars)
    if(iscellstr(vars.type))
        t=vars.type{i};
    else
        t=num2str(vars.type(i));
    end
    vars.cond{i}=[vars.cond{i} '_' t];
    vars.type{i}='mixed';
end

Stats.variables=vars;

link=Stats.probe.link;
link=link(ismember(link.type,types{1}),:);
link.type=repmat({'mixed'},height(link),1);
Stats.probe.link=link;