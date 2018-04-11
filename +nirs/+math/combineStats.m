function S = combineStats(SS)
% this function combines ChannelStats variables into a single entry by
% concatinating the fields.  



cond=unique(nirs.getStimNames(SS));
hascond = false(length(SS),length(cond));
for j=1:length(SS)
    hascond(j,:) = ismember(cond,nirs.getStimNames(SS(j)));
end

if(any(sum(hascond*1,1)>1))
    error('code only ment to combine stats with non-overlapping conditions; use nirs.modules.SubjLevelStats instead');
    return
end

S=SS(1);
S.variables=vertcat(SS.variables);
S.beta=vertcat(SS.beta);
S.covb=blkdiag(SS.covb);
S.dfe= mean(vertcat(SS.dfe));