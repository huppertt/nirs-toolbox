function a = replaceWithMissing(a,varargin)
%REPLACEWITHMISSING Insert missing data indicators into a dataset array.
%   Many functions in Statistics and Machine Learning Toolbox treat NaN as
%   representing a missing value for floating point types.  Similarly,
%   <undefined> is recognized for categorical types, the empty string for
%   cell arrays of strings, and blank strings for character arrays.
%   REPLACEWITHMISSING replaces other specified values with those
%   "standard" missing data indicators.
%
%   B = REPLACEWITHMISSING(A,...) replaces the specified numeric values,
%   categorical values, or strings in the dataset array A with NaN, <undefined>,
%   or the empty string, respectively.  Use the following parameter name/value
%   pairs to specify the values to replace:
%
%      'NumericValues'     A vector of numeric values to replace with NaN.
%      'CategoricalLevels' A string or a cell array containing strings naming
%                          the categorical levels to replace with <undefined>.
%      'Strings'           A string or a cell array containing strings to
%                          replace with the empty string.
%
%   B = REPLACEWITHMISSING(A,'DataVars',DATAVARS) replaces values only in the
%   specified dataset variables.  DATAVARS is a positive integer, a vector of
%   positive integers, a variable name, a cell array containing one or more
%   variable names, or a logical vector.  The default is all variables in A.
%
%   See also DATASET/ISMISSING, DATASET.

%   Copyright 2012-2014 The MathWorks, Inc.


pnames = {'DataVars'  'NumericValues' 'CategoricalLevels' 'Strings'};
dflts =  {       []               []                  {}        {} };
[vars,numericVals,categoricalLevels,stringVals,supplied] ...
    = dataset.parseArgs(pnames, dflts, varargin{:});

if ~supplied.DataVars
    vars = 1:a.nvars;
else
    vars = getvarindices(a,vars,false);
end
if supplied.NumericValues
    if ~isnumeric(numericVals)
        error(message('stats:dataset:replaceWithMissing:InvalidNumericValues'));
    end
end
if supplied.CategoricalLevels
    if ischar(categoricalLevels)
        categoricalLevels = cellstr(categoricalLevels);
    elseif ~iscellstr(categoricalLevels)
        error(message('stats:dataset:replaceWithMissing:InvalidCategoricalLevels'));
    end
    undef = char('');
end
if supplied.Strings
    if ischar(stringVals)
        stringVals = cellstr(stringVals);
    elseif ~iscellstr(stringVals)
        error(message('stats:dataset:replaceWithMissing:InvalidStrings'));
    end
end

a_data = a.data;
for j = 1:length(vars)
    var_j = a.data{vars(j)};
    if supplied.NumericValues && isnumeric(var_j)
        var_j(ismember(var_j,numericVals)) = NaN;
    elseif supplied.CategoricalLevels && isa(var_j,'categorical')
        var_j(ismember(var_j,categoricalLevels)) = undef;
        var_j = removecats(var_j,categoricalLevels);
    elseif supplied.Strings
        if iscellstr(var_j)
            var_j(ismember(var_j,stringVals)) = {''};
        elseif ischar(var_j)
            var_j(ismember(var_j,stringVals),:) = ' ';
        end
    end
    a_data{vars(j)} = var_j;
end
a.data = a_data;
