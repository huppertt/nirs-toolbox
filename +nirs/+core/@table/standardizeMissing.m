function a = standardizeMissing(a,missingIndicators,varargin)
%STANDARDIZEMISSING Insert standard missing data indicators into a table.
%   Many functions in MATLAB treat NaN as representing a missing value for
%   floating point types.  Similarly, <undefined> is recognized for categorical
%   types, the empty string for cell arrays of strings, and blank strings for
%   character arrays.  STANDARDIZEMISSING replaces other special values that are
%   intended to indicate missing values in a table with those "standard" missing
%   data indicators.
%
%   B = STANDARDIZEMISSING(A,INDICATORS) replaces any values specified in
%   INDICATORS that appear in floating point, categorical, or string variables
%   in the table A with NaN, <undefined>, or the empty string, respectively.
%   INDICATORS is a numeric vector, a string, or a cell array containing numeric
%   values and strings. STANDARDIZEMISSING checks floating point variables in A
%   against numeric values from INDICATORS, and string and categorical variables
%   in A against strings from INDICATORS.
%
%   Note: STANDARDIZEMISSING cannot replace values in an integer variable with NaN
%   because integer variables cannot store NaN.
%
%   B = STANDARDIZEMISSING(A,INDICATORS,'DataVariables',DATAVARS) replaces values
%   only in the specified table variables.  DATAVARS is a positive integer, a
%   vector of positive integers, a variable name, a cell array containing one or
%   more variable names, or a logical vector.  The default is all variables in A.
%
%   See also ISMISSING, TABLE.

%   Copyright 2012-2013 The MathWorks, Inc.

import matlab.internal.tableUtils.isstring

numericIndicators = [];
stringIndicators = {};
if isnumeric(missingIndicators)
    numericIndicators = missingIndicators(:);
elseif isstring(missingIndicators)
    stringIndicators = { missingIndicators }; % preserve trailing spaces
elseif iscell(missingIndicators)
    numericIndicators = cell2mat(missingIndicators(cellfun(@isnumeric,missingIndicators)));
    stringIndicators = missingIndicators(cellfun(@isstring,missingIndicators));
    if numel(numericIndicators) + numel(stringIndicators) < numel(missingIndicators)
        error(message('MATLAB:table:standardizeMissing:InvalidMissingIndicators'));
    end
else
    error(message('MATLAB:table:standardizeMissing:InvalidMissingIndicators'));
end
haveNumericIndicators = ~isempty(numericIndicators);
haveStringIndicators = ~isempty(stringIndicators);

if haveStringIndicators
    charIndicators = deblank(stringIndicators); % remove trailing spaces
    categoricalIndicators = strtrim(stringIndicators); % remove leading/trailing spaces
    % '<undefined>' is not a valid category name, don't look for it
    categoricalIndicators(strcmp(categoricalIndicators,categorical.undefLabel)) = [];
    undef = char('');
end

pnames = {'DataVariables'};
dflts =  {     1:a.nvars };
[dataVars,supplied] = matlab.internal.table.parseArgs(pnames, dflts, varargin{:});
if supplied.DataVariables
    dataVars = getVarIndices(a,dataVars);
end

% These calls to ismember will not match NaNs in numeric, or <undefined> in
% categorical. But that's OK since those elements are already standardized.

a_data = a.data;
for j = dataVars
    var_j = a_data{j};
    if haveNumericIndicators && isnumeric(var_j)
        if isfloat(var_j)
            % Only checking floats, so no chance of ismember erroring for mixed integers
            var_j(ismember(var_j,numericIndicators)) = NaN;
        end
    elseif haveStringIndicators
        if iscellstr(var_j)
            % Strict matching, no whitespace removed.
            var_j(ismember(var_j,stringIndicators)) = {''};
        elseif ischar(var_j)
            % Match without considering trailing whitespace. In particular, any empty or
            % all-blanks string matches any other.
            cvar_j = cellstr(var_j);  % remove trailing spaces
            var_j(ismember(cvar_j,charIndicators),:) = ' ';
        elseif isa(var_j,'categorical')
            % Match without considering leading/trailing whitespace.
            var_j(ismember(var_j,categoricalIndicators)) = undef;
            var_j = removecats(var_j,categoricalIndicators);
        end
    end
    a_data{j} = var_j;
end
a.data = a_data;
