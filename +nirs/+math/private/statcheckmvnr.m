function [NumSamples, NumSeries, NumParams, Y, X, goodrows] = ...
    statcheckmvnr(Y, X, Param, Covar, IgnoreNaNs)
%STATCHECKMVNR Argument checking function for mvregress.

%    Copyright 2006-2013 The MathWorks, Inc.


if nargin < 5
    IgnoreNaNs = false;
end
if nargin < 4
    Covar = [];
end
if nargin < 3
    Param = [];
end
if (nargin < 2)
    error(message('stats:statcheckmvnr:MissingInputArg'));
end

[NumSamples, NumSeries] = size(Y);

if iscell(X)
    X = X(:);
    [m, n] = size(X);

    if m == 1
        SingleX = true;
    elseif m == NumSamples
        SingleX = false;
    else
        error(message('stats:statcheckmvnr:BadXCell'));
    end

    if n > 1
        error(message('stats:statcheckmvnr:InvalidDesignArray'));
    end

    [n, NumParams] = size(X{1});

    if n ~= NumSeries
        error(message('stats:statcheckmvnr:InconsistentDims'));
    end

    if ~SingleX
        for k = 1:NumSamples
            if ~all(isequal(size(X{k}), [NumSeries, NumParams]))
                error(message('stats:statcheckmvnr:InconsistentCellDims', k, NumSeries, NumParams));
            end
        end
    end
else
    [n, NumParams] = size(X);

    if n ~= NumSamples
        error(message('stats:statcheckmvnr:InconsistentXRows', NumSamples));
    end
end

if any(sum(isnan(Y),1) == NumSamples)
    error(message('stats:statcheckmvnr:TooManyNaNs'));
end

if any(any(isinf(Y)))
    error(message('stats:statcheckmvnr:InfiniteYValue'));
end

if iscell(X)
    if SingleX
        if any(any(isnan(X{1})))
            error(message('stats:statcheckmvnr:NoMissingX'));
        end

        if any(any(isinf(X{1})))
            error(message('stats:statcheckmvnr:InfiniteValue'));
        end

        r = rank(X{1});
        if (NumSeries < NumParams) || (r ~= NumParams)
            error(message('stats:statcheckmvnr:InvalidDesignRank'));
        end
    else
        for k = 1:NumSamples
            if any(any(isinf(X{k})))
                error(message('stats:statcheckmvnr:InfiniteValueCell',k))
            end
        end
    end
else
    if any(any(isinf(X)))
        error(message('stats:statcheckmvnr:InfiniteValue'));
    end
end

if ~iscell(X) && NumSeries~=1
    if ~isempty(Param) && ~all(isequal(size(Param), [NumParams, NumSeries]))
        error(message('stats:statcheckmvnr:InconsistentBetaDims2', NumParams, NumSeries));
    end
elseif ~isempty(Param) && ~all(isequal(size(Param), [NumParams, 1]))
    error(message('stats:statcheckmvnr:InconsistentBetaDims', NumParams));
end

if ~isempty(Covar) && ~all(isequal(size(Covar), [NumSeries, NumSeries]))
    error(message('stats:statcheckmvnr:InconsistentCovarDims', NumSeries, NumSeries));
end

% Remove rows with too many NaNs
if IgnoreNaNs
    goody = ~any(isnan(Y),2);
else
    goody = ~all(isnan(Y),2);
end
if ~iscell(X)
    goodx = ~any(isnan(X),2);
elseif ~isscalar(X)
    X = X(:);
    goodx = goody;
    for j=1:numel(X)
        if any(any(isnan(X{j})))
            goodx(j) = false;
        end
    end
else
    goodx = true;  % scalar cell case was checked earlier
end

if ~any(goodx)
    error(message('stats:statcheckmvnr:InvalidDesign'));
end

goodrows = goody & goodx;
if ~all(goodrows)
    Y = Y(goodrows,:);
    NumSamples = size(Y,1);
    if ~isscalar(X)
        X = X(goodrows,:);
    end
end
