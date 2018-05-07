function Stats = conjunction(StatsIn,conditions,type,combine)
% type = {'global','Conj'}
% Global does not correct for FDR and not a proper logical AND.  See 
% Nichols et al Valid conjunction inference with the minimum statistic. NeuroImage 25 (2005) 653 ? 660 


if(nargin<3)
    type='conj';
end

if(nargin<2 || isempty(conditions))
    conditions = unique(nirs.getStimNames(StatsIn));
end

% combine flag merges HbO and HbR (or whatever types you have).
if(nargin<4)
    combine=false; 
end

StatsIn=StatsIn.sorted;

% first, lets create all the tests
T=[]; dfe=[];
for i=1:length(StatsIn)
    for j=1:length(conditions)
        if(ismember(conditions{j},StatsIn(i).conditions))
            S=StatsIn(i).ttest(conditions{j});
            T=[T S.tstat];
            dfe=[dfe S.dfe];
        end
    end
end

variab=StatsIn.variables;
if(combine)
    link=StatsIn.probe.link;
    types=unique(link.type);
    T2=[]; dfe2=[];
    for i=1:length(types)
        lst=find(ismember(link.type,types(i)));
        T2=[T2 T(lst,:)];
        dfe2=[dfe dfe];
    end
    T=T2;
    dfe=dfe2;
    lst=find(ismember(variab.type,types(1)));
    variab=variab(lst,:);
    variab.type = repmat(cellstr('joint'),height(variab),1);
end


Name=[char(8898) '{'];
for i=1:length(conditions)
    Name=[Name conditions{i} ' '];
end
Name=[Name(1:end-1) '}'];


[minT,id] = min(abs(T),[],2); 
dfe=dfe(id)';
for i=1:size(T,1)
    minT(i)=T(i,id(i));  % get the sign back
end

if(strcmp(lower(type),'global'))
     p = (2*tcdf(-abs(minT),dfe)).^size(T,2);
elseif(strcmp(lower(type),'conj'))
    p = 2*tcdf(-abs(minT),dfe);
else
    error('unknown null; must be {global, conj}');
    return;
end



tstat = minT;
q = reshape( nirs.math.fdr( p(:) )', size(p) );


variab=variab(ismember(variab.cond,variab.cond{1}),:);


Stats = [variab table(tstat, dfe, p, q)];
Stats.cond=repmat(cellstr(Name),height(Stats),1);

S=table;

flds=Stats.Properties.VariableNames;
for i=1:length(flds)
    S=[S table(Stats.(flds{i}),'VariableNames',{flds{i}})];
end

%S=setProbe(S,StatsIn(1).probe);
Stats=S;

end
