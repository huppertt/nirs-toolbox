function [tbl,lst] = find_subject_outlier(SubjStats)


t=[];
for i=1:length(SubjStats)
    subj=repmat(cellstr(['s' num2str(i)]),length(SubjStats(i).beta),1);
    t=[t; SubjStats(i).table table(subj)];
end

t=[t table(strcat(num2str(t.source),repmat('-',height(t),1), num2str(t.detector)),'VariableNames',{'srcdet'})];

lm=fitlm(t,'beta ~ subj+cond:type:srcdet', 'RobustOpts','on');

for i=1:length(SubjStats)
    rSSE(i,1)=sqrt(sum(lm.Residuals.Studentized(ismember(t.subj,cellstr(['s' num2str(i)]))).^2));
    leverage(i,1)=sum(lm.Diagnostics.Leverage(ismember(t.subj,cellstr(['s' num2str(i)]))));
    cooksD(i,1)=sqrt(sum(lm.Diagnostics.CooksDistance(ismember(t.subj,cellstr(['s' num2str(i)]))).^2));
end

p=1-tcdf(zscore(rSSE),length(SubjStats)-1);
q=nirs.math.fdr(p);

tbl = [nirs.createDemographicsTable(SubjStats) table(rSSE,cooksD,leverage,p,q)];
lst=find(tbl.q<0.05);