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

basis=SS(1).basis;
for i=2:length(SS)
    b=SS(i).basis;
    for j=1:b.base.count
        basis.base(b.base.keys{j})=b.base(b.base.keys{j});
    end
     for j=1:b.stim.count
        basis.stim(b.stim.keys{j})=b.stim(b.stim.keys{j});
    end
end


S=SS(1);
S.basis=basis;
S.variables=vertcat(SS.variables);
S.beta=vertcat(SS.beta);
S.covb=blkdiag(SS.covb);
S.dfe= mean(vertcat(SS.dfe));