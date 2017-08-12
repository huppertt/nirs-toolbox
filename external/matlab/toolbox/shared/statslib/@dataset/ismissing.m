function tf = ismissing(a,varargin)
%ISMISSING Find elements in a dataset array that contain missing values.
%   I = ISMISSING(A) returns a logical array I that indicates which elements
%   in the dataset array A contain a missing value.  B = A(~ANY(I,2),:)
%   contains only complete observations, and B = A(:,~ANY(I,1)) contains
%   only variables with no missing values.
%
%   ISMISSING recognizes NaN as indicating missing data for floating point
%   types, <undefined> for categorical types, the empty string for cell arrays
%   of strings, and blank strings for character arrays.  ISMISSING ignores
%   other variables.
%
%   Use the following parameter name/value pairs to specify additional values to
%   treat as indicating missing data:
%
%      'NumericTreatAsMissing'  A vector of numeric values to treat as missing
%                               value indicators in floating point variables.
%                               ISMISSING always treats NaN as a missing value.
%      'StringTreatAsMissing'   A string or a cell array containing strings to
%                               treat as missing value indicators in string
%                               variables.  ISMISSING always treats the empty
%                               string as a missing value.
%
%   See also DATASET/REPLACEWITHMISSING, ISNAN, ISEMPTY, ISUNDEFINED, DATASET.

%   Copyright 2012 The MathWorks, Inc.


pnames = {'NumericTreatAsMissing' 'StringTreatAsMissing'};
dflts =  {                    []                     {} };
[numericTreatAsMissing,stringTreatAsMissing,supplied] ...
    = dataset.parseArgs(pnames, dflts, varargin{:});

if supplied.NumericTreatAsMissing
    if ~isnumeric(numericTreatAsMissing)
        error(message('stats:dataset:ismissing:InvalidNumericTreatAsMissing'));
    end
end
if supplied.StringTreatAsMissing
    if ischar(stringTreatAsMissing)
        stringTreatAsMissing = cellstr(stringTreatAsMissing);
    elseif ~iscellstr(stringTreatAsMissing)
        error(message('stats:dataset:ismissing:InvalidStringTreatAsMissing'));
    end
end

tf = false(a.nobs,a.nvars);
a_data = a.data;
for j = 1:a.nvars
    var_j = a_data{j};
    tf(:,j) = statslib.internal.hasMissingVal(var_j(:,:));
    if supplied.NumericTreatAsMissing && isnumeric(var_j)
        tf(:,j) = tf(:,j) | any(ismember(var_j(:,:),numericTreatAsMissing),2);
    elseif supplied.StringTreatAsMissing && (iscellstr(var_j) || ischar(var_j))
        tf(:,j) = tf(:,j) | any(ismember(var_j(:,:),stringTreatAsMissing),2);
    end
end
