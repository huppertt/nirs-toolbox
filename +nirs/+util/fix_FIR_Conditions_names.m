function Stats=fix_FIR_Conditions_names(Stats)

if(length(Stats)>1)
    for i=1:length(Stats)
        Stats(i)=nirs.util.fix_FIR_Conditions_names(Stats(i)); 
    end
    return
end


conds=Stats.variables.cond;
condsNew=conds;
for i=1:length(conds)
   num=extractBefore(extractAfter(conds{i},':'),':');
   condsNew{i}=[extractBefore(conds{i},':') ':' ...
       extractAfter(extractAfter(conds{i},':'),':') ':' num];
end 

Stats.variables.cond=condsNew;