function terms = model2terms(modelStr,nvars,includeIntercept,treatAsCategorical)
%MODEL2TERMS Create a terms matrix from a model string
%   TERMS = MODEL2TERMS(MODELSTR,NVARS) returns a terms matrix representing
%   the model specified by MODELSTR and NVARS.  MODELSTR is a string
%   specifying the type of model.  Valid strings are 'constant', 'linear',
%   'interactions', 'purequadratic', and 'quadratic'.  NVARS is a non-negative
%   integer indicating how many variables are to be included in the model.
%   TERMS is an NTERMS-by-NVARS matrix with rows corresponding to terms in the
%   model, columns corresponding to variables in the model.  TERMS(I,J)
%   contains the power of the J-th variable in the I-th term.  A row of zeros
%   corresponds to the intercept term.
%
%   MODELSTR can also be a string of the form 'polyMN...', where [M N ...] is
%   a sequence of P numeric digits representing a P-variate polynomial in
%   which the model contains main and interaction terms for the 1st through
%   M-th powers of the first variable, 1st though N-th powers of the second
%   variable, etc.  Only those interaction terms whose powers sum to MAX([M N
%   ...}) or less are included in the model.  Variables whose specified power
%   is 0 are effectively excluded from the model.
%
%   TERMS = MODEL2TERMS(MODELSTR,WHICHVARS), where WHICHVARS is a logical
%   vector, creates an NTERMS-by-LENGTH(WHICHVARS) terms matrix TERMS.
%   WHICHVARS indicates which variables are included in the model.
%
%   TERMS = MODEL2TERMS(MODELSTR,WHICHVARS,INTERCEPT) creates the terms matrix
%   TERMS with or without an intercept term, according to the logical value
%   INTERCEPT.
%
%   TERMS = MODEL2TERMS(MODELSTR,WHICHVARS,INTERCEPT,CATVARS) treats the
%   variables indicated by the logical vector CATVARS as categorical
%   variables.  The model will include terms only for main effects and
%   interactions for these variables.
%
%   Examples:
%      >> model2terms('interactions',3)
%      ans =
%           0     0     0
%           1     0     0
%           0     1     0
%           0     0     1
%           1     1     0
%           1     0     1
%           0     1     1
%      >> model2terms('poly12',3,true,[true false true])
%      ans =
%           0     0     0
%           1     0     0
%           0     0     1
%           1     0     1
%           0     0     2

%   Copyright 2011 The MathWorks, Inc.

[c,startLoc,endLoc] = regexp(lower(modelStr),'poly(\d*)','tokens');
polyStr = ( isscalar(c) && (startLoc == 1) && (endLoc == length(modelStr)) );
    
if islogical(nvars) % model2terms(modelstr,whichvars,...)
    whichVars = nvars;
    nvars = length(whichVars);
else                % model2terms(modelstr,nvars,...)
    whichVars = true(1,nvars);
end

if nargin < 3, includeIntercept = true; end
if nargin < 4, treatAsCategorical = false(1,nvars); end

if polyStr
    powers = str2num(c{1}{1}');
    nincluded = sum(whichVars);
    if length(powers) ~= nincluded
        error(message('stats:classreg:regr:modelutils:BadLength', modelStr));
    end
    maxPower = max(powers);
    if nvars == 1
        terms = ((1-includeIntercept):maxPower)';
    else
        % Add zero powers for variables that are not included
        tmp = powers; powers = zeros(1,nvars); powers(whichVars) = tmp;
        
        % No terms above linear/interaction for categorical variables
        powers(treatAsCategorical) = max(powers(treatAsCategorical),1);
        
        % Create all terms as a cartesian product
        powers = arrayfun(@(n)0:n,powers,'UniformOutput',false);
        [powers{1:nvars}] = ndgrid(powers{:});
        powers = cellfun(@(c)c(:),powers,'UniformOutput',false);
        terms = [powers{:}];
        
        % Sort by increasing sum of powers
        [sumPowers,ord] = sort(sum(terms,2));
        terms = terms(ord,:);
        
        % Remove terms we don't want
        terms = terms(sumPowers<=maxPower,:);
        if ~includeIntercept, terms(1,:) = []; end
    end
    
else
    switch lower(modelStr)
    case 'constant',      linear = false; interactions = false; quadratic = false;
    case 'linear',        linear = true;  interactions = false; quadratic = false;
    case 'interactions',  linear = true;  interactions = true;  quadratic = false;
    case 'purequadratic', linear = true;  interactions = false; quadratic = true;
    case 'quadratic',     linear = true;  interactions = true;  quadratic = true;
    otherwise
        modelStrings = {'constant' 'linear' 'interactions' 'quadratic' 'purequadratic' 'poly*'};
        error(message('stats:classreg:regr:modelutils:BadModel', internal.stats.listStrings( modelStrings )));
    end

    icpt = zeros(includeIntercept,nvars);
    if linear
        linear = eye(nvars); linear = linear(whichVars,:);
    else
        linear = zeros(0,nvars);       
    end
    if interactions
        [rep1,rep2] = allpairs(find(whichVars));
        ninteractions = length(rep1);
        interactions = zeros(ninteractions,nvars);
        interactions(sub2ind(size(interactions),1:ninteractions,rep1)) = 1;
        interactions(sub2ind(size(interactions),1:ninteractions,rep2)) = 1;
    else
        interactions = zeros(0,nvars);
    end
    if quadratic
        whichQuadraticVars = whichVars & ~treatAsCategorical; % no quadratic for categorical
        quadratic = 2*eye(nvars); quadratic = quadratic(whichQuadraticVars,:);
    else
        quadratic = zeros(0,nvars);
    end

    terms = [icpt; linear; interactions; quadratic];
end


%-----------------------------------------------------------------------------
function [rep1,rep2] = allpairs(i,diag)
if nargin<2 || ~diag
    [r,c] = find(tril(ones(length(i)),-1));
else
    [r,c] = find(tril(ones(length(i)),0));
end
i = i(:)'; % force row vector
rep1 = i(1,c);
rep2 = i(1,r);
