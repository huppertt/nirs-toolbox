function crit = modelcriterion(stats,type,asStruct)
% MODELCRITERION Internal utility to compute criteria such as AIC

%   Copyright 2011 The MathWorks, Inc.

if nargin < 2
    type = {'AIC'};
end
[tf,type] = internal.stats.isStrings(type);
if ~tf
    error(message('stats:classreg:regr:modelutils:CriterionNotString'));
end

if any(strcmpi('all',type))
    type = {'AIC' 'AICc' 'BIC' 'CAIC'};
end

try
    L = stats.LogLikelihood;
    k = stats.NumCoefficients;
    n = stats.NumObservations;
catch me
    error(message('stats:classreg:regr:modelutils:BadStatsStructure'));
end

crit = zeros(1,length(type));
for i = 1:length(type)
    switch lower(type{i})
    case 'aic'
        crit(i) = 2*(k - L);
        type{i} = 'AIC';
    case 'aicc'
        crit(i) = 2*(k - L + k*(k+1)/(n-k-1));
        type{i} = 'AICc';
    case 'caic'
        crit(i) = k*(log(n)+1) - 2*L;
        type{i} = 'CAIC';
    case 'bic'
        crit(i) = k*log(n) - 2*L;
        type{i} = 'BIC';
    otherwise
        error(message('stats:classreg:regr:modelutils:UnknownModelCriterion',type{i}));
    end
end

if nargin > 2 && asStruct
    crit = cell2struct(num2cell(crit),type,2);
end
