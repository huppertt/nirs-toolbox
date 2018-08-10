function tf = ismissing(t,missingIndicators)
%ISMISSING Find elements in a table that contain missing values.
%   I = ISMISSING(T) returns a logical array I that indicates which elements
%   in the table T contain a missing value.  T(~ANY(I,2),:) contains only
%   complete rows, and T(:,~ANY(I,1)) contains only variables with no
%   missing values.
%
%   ISMISSING recognizes NaN as indicating missing data for floating point
%   types, <undefined> for categorical types, the empty string for cell arrays
%   of strings, and blank strings for character arrays.  ISMISSING ignores
%   other variables.
%
%   I = ISMISSING(T,INDICATORS) treats the values in INDICATORS as missing value
%   indicators.  INDICATORS is a double numeric vector, a string, or a cell array
%   containing double numeric values and strings.  ISMISSING checks numeric variables
%   in T against numeric values from INDICATORS, and string and categorical
%   variables in T against strings from INDICATORS. You must include NaN, the
%   empty string, or <undefined> in V to have ISMISSING recognize them.
%
%   Note: Integer variables cannot store NaN, therefore you must use a special
%   integer value that is otherwise unused to indicate missing integer data.
%
%   See also STANDARDIZEMISSING, ISNAN, ISEMPTY, ISUNDEFINED, TABLE.

%   Copyright 2012-2014 The MathWorks, Inc.

import matlab.internal.tableUtils.isstring
import matlab.internal.tableUtils.hasMissingVal
import matlab.internal.tableUtils.matricize

indicatorsSupplied = (nargin > 1);
numericIndicators = [];
stringIndicators = {};
% Require double to allow comparison to all numeric types, including integers.
if indicatorsSupplied
    if isa(missingIndicators,'double')
        numericIndicators = missingIndicators(:);
    elseif isstring(missingIndicators)
        stringIndicators = { missingIndicators }; % preserve trailing spaces
    elseif iscell(missingIndicators)
        numericIndicators = cell2mat(missingIndicators(cellfun(@(x)isa(x,'double'),missingIndicators)));
        stringIndicators = missingIndicators(cellfun(@isstring,missingIndicators));
        if numel(numericIndicators) + numel(stringIndicators) < numel(missingIndicators)
            error(message('MATLAB:table:ismissing:InvalidMissingIndicators'));
        end
    else
        error(message('MATLAB:table:ismissing:InvalidMissingIndicators'));
    end
end
haveNumericIndicators = ~isempty(numericIndicators);
haveStringIndicators = ~isempty(stringIndicators);
    
% ismember for numeric will not match up NaNs, so if NaN is specified in the
% list of indicators, look specifically for it in the variable.
if haveNumericIndicators
    checkForNaNs = any(isnan(numericIndicators(:)));
end
% ismember for categorical will not match '<undefined>' or empty strings to
% undefined elements in the variable, so if those are specified in the list of
% indicators, look specifically for <undefined> in the variable
if haveStringIndicators
    charIndicators = deblank(stringIndicators); % remove trailing spaces
    categoricalIndicators = strtrim(stringIndicators); % remove leading/trailing spaces
    checkForUndefineds = cellfun('isempty',categoricalIndicators) ...
                            | any(strcmp(categorical.undefLabel,categoricalIndicators(:)));
end

tf = false(t.nrows,t.nvars);
t_data = t.data;
for j = 1:t.nvars
    var_j = matricize(t_data{j});
    if ~indicatorsSupplied
        tf(:,j) = hasMissingVal(var_j);
    elseif isnumeric(var_j) && haveNumericIndicators
        if checkForNaNs
            tf(:,j) = any(ismember(var_j,numericIndicators),2) | any(isnan(var_j),2);
        else
            tf(:,j) = any(ismember(var_j,numericIndicators),2);
        end            
    elseif haveStringIndicators
        if iscellstr(var_j)
            % Strict matching, no whitespace removed. In particular, empty or all-blanks
            % strings must match exactly.
            tf(:,j) = any(ismember(var_j,stringIndicators),2);
        elseif ischar(var_j)
            % Match without considering trailing whitespace. In particular, any empty or
            % all-blanks string matches any other.
            var_j = matricize(cellstr(t_data{j})); % removes trailing spaces
            tf(:,j) = any(ismember(var_j,charIndicators),2);
        elseif isa(var_j,'categorical')
            % Match without considering leading/trailing whitespace.
            if checkForUndefineds
                tf(:,j) = any(ismember(var_j,categoricalIndicators),2) | any(isundefined(var_j),2);
            else
                tf(:,j) = any(ismember(var_j,categoricalIndicators),2);
            end
        end
    end
end
