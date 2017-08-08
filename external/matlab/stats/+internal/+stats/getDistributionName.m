function name = getDistributionName(name)
%getDistributionName Get standard distribution name from shortened name.
%   NAME = internal.stats.getDistributionName(NAME) take a text string
%   representing a distribution name, possibly abbreviated, and returns it
%   in a standard form. Unrecognized names get returned unchanged.

%   Copyright 2014 The MathWorks, Inc.

distNames = {'beta', 'binomial', 'chi-square', 'extreme value', ...
    'exponential', 'f', 'gamma', 'generalized extreme value', ...
    'generalized pareto', 'geometric', 'hypergeometric', ...
    'lognormal', 'negative binomial', 'noncentral f', ...
    'noncentral t', 'noncentral chi-square', 'normal', 'poisson', ...
    'rayleigh', 't', 'discrete uniform', 'uniform', 'weibull'};

i = find(strncmpi(name, distNames, length(name)));
if numel(i) > 1
    throwAsCaller(MException('stats:random:AmbiguousDistribution',name));
elseif numel(i) == 1
    name = distNames{i};
    return
end

% it may be an abbreviation that doesn't partially match the name
name = lower(name);

switch name
    case 'ev'
        name = 'extreme value';
    case 'gev'
        name = 'generalized extreme value';
    case 'gp'
        name = 'generalized pareto';
    case 'nbin'
        name = 'negative binomial';
    case 'wbl'
        name = 'weibull';
    case {'chi2','chisquare'}
        name = 'chi-square';
    case 'hyge'
        name = 'hypergeometric';
    case 'ncf'
        name = 'noncentral f';
    case 'nct'
        name = 'noncentral t';
    case 'ncx2'
        name = 'noncentral chi-square';
    case 'unid'
        name = 'discrete uniform';
end