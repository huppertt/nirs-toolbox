function r = rsquared(stats,type,asStruct)
% RSQUARED Internal utility for displaying r-squared values.

%   Copyright 2011-2013 The MathWorks, Inc.

if nargin < 2
    type = {'Ordinary'};
end
[tf,type] = internal.stats.isStrings(type);
if ~tf
    error(message('stats:classreg:regr:modelutils:InvalidRsquaredType'));
end

try
    SSE = stats.SSE;
    SST = stats.SST;
catch me
    error(message('stats:classreg:regr:modelutils:InvalidRsquaredStruct'));
end

if any(strcmpi('all',type))
    type = {'Ordinary' 'Adjusted' 'LLR' 'Deviance' 'AdjGeneralized'};
end

r = zeros(1,length(type));
for i = 1:length(type)
    switch lower(type{i})
    case 'ordinary'
        r(i) = 1 - SSE./SST;
        type{i} = 'Ordinary';
    case 'adjusted'
        try
            DFE = stats.DFE;
            DFT = stats.NumObservations - 1;
        catch me
            error(message('stats:classreg:regr:modelutils:InvalidRsquaredStructAdj'));
        end
        r(i) = 1 - (SSE./SST)*(DFT./DFE);
        type{i} = 'Adjusted';
    case 'llr'
        try
            L = stats.LogLikelihood;
            L0 = stats.LogLikelihoodNull;
        catch me
            error(message('stats:classreg:regr:modelutils:InvalidRsquaredStructLLR'));
        end
        r(i) = 1 - L./L0;
        type{i} = 'LLR';
    case 'deviance'
        try
            D = stats.Deviance;
            D0 = stats.DevianceNull;
        catch me
            error(message('stats:classreg:regr:modelutils:InvalidRsquaredStructDev'));
        end
        r(i) = (D0-D)./D0;
        type{i} = 'Deviance';
    case 'adjgeneralized'
        try
            L = stats.LogLikelihood;
            L0 = stats.LogLikelihoodNull;
            N = stats.NumObservations;
        catch me
            error(message('stats:classreg:regr:modelutils:InvalidRsquaredStructAGR'));
        end
        r(i) = (1 - exp(2*(L0 - L)/N)) ./ (1 - exp(2*L0/N));
        type{i} = 'AdjGeneralized';
    otherwise
        error(message('stats:classreg:regr:modelutils:UnknownRsquaredType',type{i}));
    end
end

if nargin > 2 && asStruct
    r = cell2struct(num2cell(r),type,2);
end
