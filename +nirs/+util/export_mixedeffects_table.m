function varargout=export_mixedeffects_table(S,centerVars)
%This function does the same code used in the mixed effects models (nirs.modules.MixedEffects), but is intended
%to export the variables as a table that can be used in SPSS/R/etc

if(nargin<2)
    centerVars=true;
end

formula='beta~-1+cond';

% demographics info
demo = nirs.createDemographicsTable( S );
            
% center numeric variables
if centerVars
    n = demo.Properties.VariableNames;
    for i = 1:length(n)
        if all( isnumeric( demo.(n{i}) ) )
            demo.(n{i}) = demo.(n{i}) - nanmean( demo.(n{i}) );
        end
    end
end
            
% preallocate group stats
G = nirs.core.ChannelStats();
%% loop through files
W = sparse([]);
iW = sparse([]);

b = [];

LstV=[];
vars = table();
for i = 1:length(S)

    lstValid=~isnan(S(i).tstat);
    LstV=[LstV; lstValid];
    % coefs
    if ~isempty(strfind(formula(1:strfind(formula,'~')-1),'tstat'))
                    b = [b; S(i).tstat];
                else
                    b = [b; S(i).beta];
                end
    
    % whitening transform
    L = chol(S(i).covb,'upper');
    W = blkdiag(W,pinv(L));

    % table of variables
    file_idx = repmat(i, [height(S(i).variables) 1]);

    if(~isempty(demo))
        vars = [vars;
            [table(file_idx) S(i).variables repmat(demo(i,:), [height(S(i).variables) 1])]
            ];
    else
        vars = [vars; ...
            [table(file_idx) S(i).variables]];
    end
end

% sort
if(~ismember('source',vars.Properties.VariableNames) & ...
        ismember('ROI',vars.Properties.VariableNames))
    [vars, idx] = nirs.util.sortrows(vars, {'ROI', 'type'});

    % list for first source
    [sd, ~,lst] = nirs.util.uniquerows(table(vars.ROI, vars.type));
    sd.Properties.VariableNames = {'ROI', 'type'};

elseif(ismember('NameKernel',vars.Properties.VariableNames))
    [vars, idx] = nirs.util.sortrows(vars, {'NameKernel', 'type'});

    % list for first source
    [sd, ~,lst] = nirs.util.uniquerows(table(vars.NameKernel, vars.type));
    sd.Properties.VariableNames = {'NameKernel', 'type'};

else

    [vars, idx] = nirs.util.sortrows(vars, {'source', 'detector', 'type'});

    % list for first source
    [sd, ~,lst] = nirs.util.uniquerows(table(vars.source, vars.detector, vars.type));
    sd.Properties.VariableNames = {'source', 'detector', 'type'};
end


%% design mats
for c = 1:height(vars)
    block_ind = strfind(vars.cond{c},'◄');
    if ~isempty(block_ind)
        vars.cond{c} = vars.cond{c}(1:block_ind-2);
    end
end

tmp = vars(lst == 1, :);
beta = randn(size(tmp,1), 1);

nRE=max(1,length(strfind(formula,'|')));
warning('off','stats:LinearMixedModel:IgnoreCovariancePattern');

conds=unique(tmp.cond);
for i=1:length(conds);
    tmp.(conds{i})=1*ismember(tmp.cond,conds{i});
end;

NoInter=[];
if(~isempty(strfind(formula,'{')))
    % the formula has variables of no interest
    lstt=sort([strfind(formula,'{') strfind(formula,'}')]);
    cnt=1;
    for ii=1:2:length(lstt)
        NoInter{cnt}=formula([lstt(ii)+1:lstt(ii+1)-1]);
        cnt=cnt+1;
    end
    formula=strrep(formula,'{',' ');
    formula=strrep(formula,'}',' ');
end


formula=nirs.util.verify_formula([table(beta) tmp], formula,true);
respvar = formula(1:strfind(formula,'~')-1);


ww=full(diag(W));
data_tbl = [table(b(idx),ww(idx),'VariableNames',{respvar,'weights'}) vars];
varNames = data_tbl.Properties.VariableNames;

for i = 1:numel(varNames)
    var = data_tbl.(varNames{i});
    if ~isnumeric(var) && ~islogical(var)
        % Convert to categorical if not already
        if ~iscategorical(var)
            var = categorical(var);
        end
        % Reorder categories alphabetically
        data_tbl.(varNames{i}) = reordercats(var);
    end
end

data_tbl.measurement=categorical(strcat('Src',num2str(data_tbl.source),'_Det',num2str(data_tbl.detector),'_',cellstr(data_tbl.type)));

varargout{1}=data_tbl;

if(nargout==2)
    varargout{2}=W;
end