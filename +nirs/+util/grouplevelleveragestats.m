function tbl=grouplevelleveragestats(S,formula,centerVars)
% This function computes the leverage (across channels and across subjects)
% for a second level (group) model

if(length(S)==1)
    tbl=[];
    warning('must have more than one file for outlier definitions');
    return;
end

if(nargin<2 || isempty(formula))
    formula='beta ~ -1 + cond';
end

if(nargin<3)
    centerVars=true;
end

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
            
%% loop through files
W = sparse([]);
iW = sparse([]);

b = [];
vars = table();
vars_short = table();
for i = 1:length(S)
    % coefs
    if ~isempty(strfind(formula(1:strfind(formula,'~')-1),'tstat'))
        b = [b; S(i).tstat];
    else
        b = [b; S(i).beta];
    end
    
    % whitening transform
    
    [u, s, ~] = svd(S(i).covb, 'econ');
    %W = blkdiag(W, diag(1./diag(sqrt(s))) * u');
    W = blkdiag(W, pinv(s).^.5 * u');
    
    iW = blkdiag(iW, u*sqrt(s) );
    
    
    %                L = chol(S(i).covb,'upper');
    %                W = blkdiag(W,pinv(L));
    
    % table of variables
    file_idx = repmat(i, [size(S(i).beta,1) 1]);
    
    if(~isempty(demo))
        vars_short = [vars_short; ...
            [table(file_idx) S(i).variables]];
        vars = [vars;
            [table(file_idx) S(i).variables repmat(demo(i,:), [size(S(i).beta,1) 1])]
            ];
        
    else
        vars = [vars; ...
            [table(file_idx) S(i).variables]];
        vars_short=vars;
    end
end

% sort
[vars, idx] = nirs.util.sortrows(vars, {'source', 'detector', 'type'});
[vars_short] = nirs.util.sortrows(vars_short, {'source', 'detector', 'type'});


% list for first source
[sd, ~,lst] = nirs.util.uniquerows(table(vars.source, vars.detector, vars.type));
sd.Properties.VariableNames = {'source', 'detector', 'type'};

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

formula=nirs.util.verify_formula([table(beta) tmp], formula,true);
respvar = formula(1:strfind(formula,'~')-1);

lm1 = fitlme([table(beta,'VariableNames',{respvar}) tmp],formula, 'dummyVarCoding',...
    'full', 'FitMethod', 'ML', 'CovariancePattern', repmat({'Isotropic'},nRE,1));

X = lm1.designMatrix('Fixed');
Z = lm1.designMatrix('Random');

if(~isempty(Z))
    warning('random effects are ignored when computing leverage');
end

nchan = max(lst);

X = kron(speye(nchan), X);

%% put them back in the original order
vars=vars_short;
vars(idx,:) = vars;
X(idx, :)   = X;
beta        = b; % already in correct order


R = qr(X,0);
E = X/R;
h = sum((E.*E)')';

WX=W*X;
WX=sparse(WX.*(abs(WX)>eps(1)));

R = qr(WX,0);
E = (WX)/R;
hw = sum((E.*E)')';

tb=vars;
tb.source=[];
tb.detector=[];
if any(strcmp(tb.Properties.VariableNames,'ROI'))
    tb.ROI=[];
end

[utb,~,j]=unique(tb);
lev_condition=zeros(height(tb),1);
pval_condition=zeros(height(tb),1);
for i=1:height(utb)
    lev_condition(j==i)=mean(hw(j==i));
end
pval_condition=2*tcdf(-abs(zscore(lev_condition)),height(utb));




tb.type=[];
tb.cond=[];

[utb,~,j]=unique(tb);
lev_subject=zeros(height(tb),1);
pval_subject=zeros(height(tb),1);
for i=1:height(utb)
    lev_subject(j==i)=mean(hw(j==i));
    
end
pval_subject=2*tcdf(-abs(zscore(lev_subject)),height(utb));


lev_unweighted=full(h);
pval_unweighted=1-2*tcdf(-abs(zscore(lev_unweighted)),height(vars));
lev_weighted=full(hw);
pval_weighted=2*tcdf(-abs(zscore(lev_weighted)),height(vars));


tbl=[vars table(lev_unweighted,pval_unweighted,lev_weighted,pval_weighted,lev_condition,pval_condition,lev_subject,pval_subject)];


            