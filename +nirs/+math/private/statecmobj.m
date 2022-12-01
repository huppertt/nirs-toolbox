function Objective = statecmobj(Design, Data, Param, Covar, Resid, CholCovar, isdiagonal)
%STATECMOBJ Objective function for ecm multivariate regression.
%   OBJ = STATECMOBJ(X,Y,BETA,SIGMA,RESID,CHOLSIGMA,ISDIAGONAL) computes the
%   objective function for ecm (mle or covariance weighted) multivariate
%   regression for the specified parameter values.
%
%   This function requires Covar to be positive-definite.
%
% See also MVREGRESS.

%    Copyright 2006-2007 The MathWorks, Inc.


if nargin < 3
    error(message('stats:mvregresslike:MissingInputArg'));
end

Param = Param(:);

[~, NumSeries] = size(Data);

if iscell(Design) && (numel(Design) == 1)
    SingleDesign = true;
else
    SingleDesign = false;
end

havechol = nargin>=6 && ~isempty(CholCovar);
if ~havechol
    if (nargin>=7) && isdiagonal        % use faster formula for diagonal covariance
        d = diag(Covar);
        CholState = sum(d<=0);
        CholCovar = diag(sqrt(d));
    else
        [CholCovar, CholState] = chol(Covar);
    end
elseif rank(CholCovar)<size(Covar,1)
    CholState = 1;
else
    CholState = 0;
end
if CholState > 0
    error(message('stats:mvregresslike:NonPosDefCov'));
end

LogTwoPi = log(2.0 * pi);
LogDetCovar = 2.0 * sum(log(diag(CholCovar)));

Objective = 0.0;

nans = isnan(Data);
fullrows = ~any(nans,2);
partialrows = find(~fullrows)';

% Do full rows first
if nargin<5 || isempty(Resid)
    r = find(fullrows);
    n = length(r);
    if ~iscell(Design)
        Means = Design(r,:) * Param;
    elseif SingleDesign
        Means = repmat((Design{1} * Param)', n, 1);
    else
        Means = zeros(n,size(Data,2));
        for k=1:n
            Means(k,:) = (Design{r(k)} * Param)';
        end
    end
    res = Data(fullrows,:) - Means;
else
    res = Resid(fullrows,:);
end

if any(fullrows)
    Objective = - 0.5 * sum(sum((res / CholCovar).^2,2)) ...
                - 0.5 * sum(fullrows) * (LogDetCovar  + NumSeries*LogTwoPi);
end

for k = partialrows

    P = ~nans(k,:);
    Available = sum(P);
    if Available == 0
        continue;
    end

    Objective = Objective - 0.5 * Available * LogTwoPi;

    if ~iscell(Design)
        Mean = Design(k,:) * Param;
    elseif SingleDesign
        Mean = Design{1} * Param;
    else
        Mean = Design{k} * Param;
    end

    if isdiagonal && ~havechol
        SubChol = CholCovar(P,P);
    else
        [SubChol, CholState] = chol(Covar(P,P));
        if CholState > 0
            error(message('stats:mvregresslike:NonPosDefSubCov'));
        end
    end

    SubResid = SubChol' \ (Data(k,P)' - Mean(P));

    Objective = Objective - 0.5 * (SubResid' * SubResid);
    Objective = Objective - sum(log(diag(SubChol)));
end
