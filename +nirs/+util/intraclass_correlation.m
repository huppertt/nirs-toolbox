function ICCtbl=intraclass_correlation(SubjStats,demo_variables,fixed_effects)
% This function computes intraclass correlation using a hierachical 
% linear mixed effects model.
% 
% see 
% Trial and error: a hierarchical modeling approach to test-retest assessment
% Gang Chen, Daniel S. Pine, Melissa A. Brotman, Ashley R. Smith, Robert W. Cox, Simone P. Haller
% bioRxiv 2021.01.04.425305; doi: https://doi.org/10.1101/2021.01.04.425305
% Now published in NeuroImage doi: 10.1016/j.neuroimage.2021.118647
%
% https://www.biorxiv.org/content/10.1101/2021.01.04.425305v3.full

% Inputs: 
%     SubjStats- the subject level stats model from the GLM
%     demo_variables- a cell array of the names of the demographics to include in the 
%                     model.  Usually this is something like "subject" 

% Note- if you add more demo_variables the ICC will change because you are
% controlling for the variability of more things.  

if(nargin<3)
    fixed_effects={};
end
if(~iscell(fixed_effects))
    fixed_effects={fixed_effects};
end
if(~iscell(demo_variables))
    demo_variables={demo_variables};
end

tbl=SubjStats.table;

lstVar=find(ismember(tbl.Properties.VariableNames,{'source','detector'}));
[~,~,lst]=unique(tbl(:,lstVar));

% make a new categorical variable for measurements
for i=1:length(lst); 
    ml{i,1}=['meas-' num2str(lst(i))]; 
end;
tbl.measurement=ml;

se=tbl.se;

formula='beta ~ -1 + measurement:cond:type';
fixed='-1 + cond:type';
if(~isempty(fixed_effects))
    for i=1:length(fixed_effects)
        formula=[formula ':' fixed_effects{i}];
        fixed=[fixed ':' fixed_effects{i}];
    end
end

for i=1:length(demo_variables)
    formula=[formula ' + (' fixed '|' demo_variables{i} ')'];
end

model=fitlme(tbl,formula,...
    'FitMethod','REML','DummyVarCoding','full','Weights',1./se,...
    'CovariancePattern',repmat({'Diagonal'},1,length(strfind(formula,'|'))));

% The covariance info (specifically the order of the names) is stored as
% a private field in the LinearModel.  This is the easiest way to get this 
% info out without having to edit the internal matlab functions

s=evalc('disp(model);');
s=s(strfind(s,'Random effects covariance parameters (95% CIs):')+57:end-2);
ss=strsplit(s,char(10));

ICCtbl=struct;
cnt=1;
MSE=NaN;
for i=1:length(ss)
    if(contains(ss{i},'Res Std'))
         a=textscan(ss{i},'%s%s%f%f%f');
         MSE=a{3};
    end
end
for i=1:length(ss)
    if(contains(ss{i},'std'))
        for ii=20:-1:2; ss{i}=strrep(ss{i},repmat(' ',1,ii),' '); end;
        ss{i}=strrep(ss{i},'{ ','{');
        ss{i}=strrep(ss{i},' }','}');
        ss{i}=strrep(ss{i},char(39),'');
        a=textscan(ss{i},'%s%s%s%f%f%f');
        ICCtbl.Variable(cnt,1)=string(strrep(strrep(strrep(a{1},'{',''),'}',''),'''',''));
        ICCtbl.Sigma2(cnt,1)=a{4}^2;
        ICCtbl.ICC(cnt,1)=a{4}^2/(a{4}^2+MSE);
        cnt=cnt+1;
    end
end
ICCtbl.Variable(cnt,1)="Residual Error";
ICCtbl.Sigma2(cnt,1)=MSE;
ICCtbl.ICC(cnt,1)=NaN;

ICCtbl=struct2table(ICCtbl);