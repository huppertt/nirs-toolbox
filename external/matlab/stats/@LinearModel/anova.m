function tbl = anova(model,anovatype,sstype)
%ANOVA Analysis of variance
%   TBL = ANOVA(LM) displays an analysis of variance (anova) table for the
%   LinearModel LM. The table displays a 'components' anova with sums of
%   squares and F tests attributable to each term in the model, except the
%   constant term.
%
%   TBL = ANOVA(LM,'summary') displays a summary anova table with an F test
%   for the model as a whole. If there are both linear and higher-order
%   terms, there is also an F test for the higher-order terms as a group.
%   If there are replications (multiple observations sharing the same
%   predictor values), there is also an F test for lack-of-fit computed by
%   decomposing the residual sum of squares into a sum of squares for the
%   replicated observations and the remaining sum of squares.
%
%   TBL = ANOVA(LM,'components',SSTYPE) computes the specified type of sums
%   of squares. Choices 1-3 give the usual type 1, type 2, or type 3 sums
%   of squares.  The value 'h' produces sums of squares for a hierarchical
%   model that is similar to type 2, but with continuous as well as
%   categorical factors used to determine the hierarchy of terms. The
%   default is 'h'.
%
%   Example:
%       % Look at a components anova showing all terms
%       load carsmall
%       d = dataset(MPG,Weight);
%       d.Year = ordinal(Model_Year);
%       lm = LinearModel.fit(d,'MPG ~ Year + Weight + Weight^2')
%       anova(lm)
%
%       % Look at a summary anova showing groups of terms. The nonlinear
%       % group consists of just the Weight^2 term, so it has the same
%       % p-value as that term in the previous table. The F statistic
%       % comparing the residual sum of squares to a "pure error" estimate
%       % from replicated observations shows no evidence of lack of fit.
%       anova(lm,'summary')
%
%   See also LinearModel.

%   Copyright 2011-2013 The MathWorks, Inc.

if nargin < 2
    anovatype = 'components';
else
    anovatype = internal.stats.getParamVal(anovatype,...
        {'summary' 'components'},'second');
end
if nargin<3
    sstype = 'h';
end
switch(lower(anovatype))
case 'components'
    tbl = componentanova(model,sstype);
case 'summary'
    tbl = summaryanova(model);
otherwise
    error(message('stats:LinearModel:BadAnovaType'));
end


% ----------------------------------------------------------------------
function tbl = componentanova(model,sstype)


% Initialize variables for components anova
nterms = length(model.Formula.TermNames);

if sstype==3
    % Re-do model using the dummy variable coding that we require.
    oldmodel = model;
    formula = sprintf('%s ~ %s',oldmodel.Formula.ResponseName,oldmodel.Formula.LinearPredictor);
    model = LinearModel.fit(oldmodel.Variables,formula,'DummyVarCoding','effects',...
        'Robust',oldmodel.Robust,'Categorical',oldmodel.VariableInfo.IsCategorical,...
        'Exclude',~oldmodel.ObservationInfo.Subset,'Weight',oldmodel.ObservationInfo.Weights);
end

% The SS for a term is the difference in RSS between a model with the term
% and a model without the term. The specific models depend on sstype.
allmodels = getComparisonModels(model,sstype);
[uniquemodels,~,uidx] = unique(allmodels,'rows');

% Find SS and DF for each model
allSS = zeros(nterms,2);
allDF = zeros(nterms,2);
dfx = model.NumObservations - model.DFE;
for j=1:size(uniquemodels,1)
    [ssModel,dfModel] = getTermSS(model,~uniquemodels(j,:),dfx);
    t = (uidx==j);
    allSS(t) = ssModel;
    allDF(t) = dfModel;
end

ss = [diff(allSS,1,2); model.SSE];
df = [diff(allDF,1,2); model.DFE];

% Remove constant, if any
constTerm = all(model.Formula.Terms==0,2);
ss(constTerm,:) = [];
df(constTerm,:) = [];

% Compute remaining columns of the table
ms = ss ./ df;
invalidrows = (1:length(ms))' == length(ms);
f = ms./ms(end);
pval = fcdf(1./f,df(end),df); % 1-fcdf(f,df,dfe)

% Assemble table
tbl = table(ss, df, ms, ...
    internal.stats.DoubleTableColumn(f,invalidrows), ...
    internal.stats.DoubleTableColumn(pval,invalidrows), ...
    'VariableNames',{'SumSq' 'DF' 'MeanSq' 'F' 'pValue'}, ...
    'RowNames',[model.Formula.TermNames(~constTerm);'Error']);


% ----------------------------------------------------------------------
function tbl = summaryanova(model)
% Create summary anova table with rows
%
% 1   Total SS
% 2   Model SS
% 3      Linear SS        <--- These two rows are present only if
% 4      Nonlinear SS          the model has interaction or power terms
% 5   Residual SS
% 6      Lack-of-fit SS   <--- These two rows are present only if
% 7      Pure error SS         the model has replicates

ss = zeros(7,1);
df = zeros(7,1);
keep = true(7,1);

termorder = sum(model.Formula.Terms,2);

% Get information always required
hasconst = any(termorder==0);
ss(1) = model.SST;
df(1) = model.NumObservations - hasconst;

ss(2) = model.SSR;
dfx = df(1) - model.DFE;
df(2) = dfx;

ss(5) = model.SSE;
df(5) = model.DFE;

% If there are nonlinear terms, we can break ssr into pieces
if sum(termorder>1)>0 && sum(termorder==1)>0
    terminfo = getTermInfo(model);
    nonlincols = ismember(terminfo.designTerms, find(termorder>1));
    [ss(4),df(4)] = getTermSS(model,nonlincols,dfx+hasconst);
    ss(3) = ss(2) - ss(4);
    df(3) = df(2) - df(4);
else
    keep(3:4) = false;
end

% If there are replicates, we can break sse into pieces
subset = model.ObservationInfo.Subset;
[sx,ix] = sortrows(model.Design(subset,:)); % better to sort X itself?
isrep = [all(diff(sx)==0,2); false];
sx = []; % no longer needed, save space
if any(isrep)
    sspe = 0;   % sum of squares for pure error
    dfpe = 0;   % degrees of freedom for pure error
    first = 1;
    n = length(isrep);
    r = model.Residuals.Raw(subset);
    w = model.ObservationInfo.Weights(subset);
    while(first<n)
        % find stretch of replicated observations
        if ~isrep(first)
            first = first+1;
            continue;
        end
        for k=first+1:n
            if ~isrep(k)
                last = k;
                break
            end
        end
        
        % Compute pure error SS and DF contributions from this stretch
        t = ix(first:last);
        r1 = r(t);
        w1 = w(t);
        m = sum(w1.*r1) / sum(w1);
        sspe = sspe + sum(w1.*(r1-m).^2);
        dfpe = dfpe + (last-first);
        
        % Continue beyond this stretch
        first = last+1;
    end
    if dfpe==model.DFE
        % Replications but no lack-of-fit, so treat as if no reps
        isrep = false;
    end
end
if any(isrep)
    ss(7) = sspe;
    df(7) = dfpe;
    ss(6) = ss(5) - ss(7);
    df(6) = df(5) - df(7);
else
    keep(6:7) = false;
end

% Compute MS for each term
ms = ss ./ df;

% Define error terms for each F test
invalidrows = [true false false false true false true]';
mse = [NaN ms(5) ms(5) ms(5) NaN ms(7) NaN]';
dfe = [NaN df(5) df(5) df(5) NaN df(7) NaN]';

f = ms./mse;
pval = fcdf(1./f,dfe,df); % 1-fcdf(f,df,dfe)
tbl = table(ss, df, ms, ...
    internal.stats.DoubleTableColumn(f,invalidrows), ...
    internal.stats.DoubleTableColumn(pval,invalidrows), ...
    'VariableNames',{'SumSq' 'DF' 'MeanSq' 'F' 'pValue'}, ...
    'RowNames',{'Total' 'Model' '. Linear' '. Nonlinear' ...
    'Residual' '. Lack of fit' '. Pure error'});
tbl = tbl(keep,:);

obsnames = {'Total' 'Model' '. Linear' '. Nonlinear' ...
    'Residual' '. Lack of fit' '. Pure error'};
tbl = table(ss(keep), df(keep), ms(keep), ...
    internal.stats.DoubleTableColumn(f(keep),invalidrows(keep)), ...
    internal.stats.DoubleTableColumn(pval(keep),invalidrows(keep)), ...
    'VariableNames',{'SumSq' 'DF' 'MeanSq' 'F' 'pValue'}, ...
    'RowNames',obsnames(keep));



% ----------------------------------------------------------------------
function allmodels = getComparisonModels(model,sstype)
terminfo = getTermInfo(model);
termcols = terminfo.designTerms;
nterms = max(termcols);
terms = model.Formula.Terms;
allmodels = false(2*nterms,length(termcols));
continuous = ~terminfo.isCatVar;
switch(sstype)
    case 1
        % Type 1 or sequential sums of squares
        for j = 1:nterms
            allmodels(j,:)        = termcols<=j; % with term j
            allmodels(j+nterms,:) = termcols<j;  % without term j
        end
    case 3
        for j = 1:nterms
            allmodels(j,:)        = true;        % with term j
            allmodels(j+nterms,:) = termcols~=j; % without term j
        end
    case {2,'h'}
        % Strict hierarchical sums of squares
        for j = 1:nterms
            % Get vars in this term
            varsin = terms(j,:);
            
            % Take out terms higher than this one
            out = all(bsxfun(@ge,terms(:,varsin>0),terms(j,varsin>0)),2);
            
            % But for type 2, only if they match on all continuous vars
            if sstype==2
                out = out & all(bsxfun(@eq,terms(:,continuous),terms(j,continuous)),2);
            end
            t = ismember(termcols,find(~out));
            allmodels(j,:)        = t | termcols==j; % with term j
            allmodels(j+nterms,:) = t;               % without term j
        end
end

% ----------------------------------------------------------------------
function [ss,df] = getTermSS(model,termcols,dfbefore)
% Compute SS and DF for term or terms in specific cols
R = model.R;

if size(R,1)==dfbefore && model.MSE>0 
    % This is a full-rank model, so we can computes sums of squares using
    % the coefficient estimates and their covariance.
    b = model.Coefs(termcols);
    b = b(:);   % need column vector even if empty
    V = model.CoefficientCovariance(termcols,termcols);
    ss = b'*(V\b) * model.MSE;
    df = length(b);
else
    % This is a reduced-rank model with some coefficients artificially set
    % to zero and their variances set to 0. We have to operate on the R
    % matrix to compute the degrees of freedom and sum of squares.
    % If X is n-by-p, this reduces the problem to p-by-p.
    Qy = model.Qy;
    Rtol = model.Rtol;
    
    % Refactor R so q spans remaining space. df is the number of independent
    % columns. Error SS is the SS from the orthogonal part
    [q,r,~] = qr(R(:,~termcols));
    if isvector(r)
        dfafter = any(r);
    else
        dfafter = sum(abs(diag(r)) > Rtol);
    end
    ss = norm(Qy'*q(:,dfafter+1:end))^2;
    df = dfbefore - dfafter;
end
