function OS = createNonConstantDefaults(OS)
%

%CREATENONCONSTANTDEFAULTS Create the NonConstantDefaults structure
%
%   OS = CREATENONCONSTANTDEFAULTS(OS) creates the non-constant defaults
%   structures, OS.NonConstantDefaultFields, OS.NonConstantDefaults and
%   OS.IsConstantDefault.
% 
%   OS.NonConstantDefaultFields is a cell-string containing names of the
%   options where the set of default values over all the algorithms
%   contains more than one value. We define this as an option with a
%   "non-constant" default.
%   
%   OS.NonConstantDefaults contains the default values for those options
%   whose default differs between the algorithms.
%
%   OS.IsConstantDefault is a structure containing an entry for each solver
%   option. The value of each field is a boolean, which determines whether
%   that option has a constant default or not. For example, for fmincon,
%   OS.IsConstantDefault.TolX = false indicates that TolX does not have a
%   constant default for fmincon.
%
%   See the architecture specification for more information on non-constant
%   defaults.

%   Copyright 2012 The MathWorks, Inc.

% Number of algorithms
numAlgs = length(OS.AlgorithmDefaults);

% Generate the non constant default field names as a cell array

% Get a list of all the solver options
fnames = {};
for i = 1:numAlgs
    fnames = [fnames; fieldnames(OS.AlgorithmDefaults{i})]; %#ok
end
fnames = unique(fnames);

% Create a 1-by-numFnames logical array which indicates whether option
% fnames{i} has a non-constant default.
numFnames = length(fnames);
isNCD = false(1, numFnames);

% Loop through each option. Check the default values for each algorithm for
% the given option. If more than one default value exists for the option,
% then mark the option as having a non-constant default.
for i = 1:numFnames
    firstVal = true;
    for j = 1:numAlgs
        if isfield(OS.AlgorithmDefaults{j}, fnames{i})
            thisVal = OS.AlgorithmDefaults{j}.(fnames{i});
            if firstVal
                refVal = thisVal;
                firstVal = false;
            else
                % Check to see if thisVal is equal to refVal. We are
                % testing default values of the options in the Optimization
                % solvers here. Currently default option values take one of
                % the following types
                % * String
                % * Numeric
                % * Named Function handle (e.g @gacreationuniform)
                if ischar(refVal)
                    % Equality test for strings
                    refValEqualToThisVal = strcmp(refVal, thisVal);
                else
                    % Equality test for numeric values and named function
                    % handles. We also use this to catch any other data
                    % type for now.
                    refValEqualToThisVal = isequal(refVal, thisVal);
                end
                    
                % If a second default value exists for the option, mark the
                % option as having a non-constant default and move to the
                % next option.
                if ~refValEqualToThisVal
                    isNCD(i) = true;
                    break
                end
            end
        end
    end
end
OS.NonConstantDefaultFields = fnames(isNCD);

% Create the IsConstantDefault structure
for i = 1:length(fnames)
    OS.IsConstantDefault.(fnames{i}) = true;
end
for i = 1:length(OS.NonConstantDefaultFields)
    OS.IsConstantDefault.(OS.NonConstantDefaultFields{i}) = false;
end

% Create the NonConstantDefaults structure
for i = 1:length(OS.NonConstantDefaultFields)
    thisDefaults = cell(1, numAlgs);
    idxUndefined = false(1, numAlgs);
    
    % Generate the NonConstantDefaults cell for the current field
    for j = 1:numAlgs
        if isfield(OS.AlgorithmDefaults{j}, OS.NonConstantDefaultFields{i})
            thisDefaults{j} = OS.AlgorithmDefaults{j}.(OS.NonConstantDefaultFields{i});
        else
            idxUndefined(j) = true;
        end
    end
    
    % Insert the correct undefined value
    if any(idxUndefined)
       idxDefined = find(~idxUndefined);
       thisVal = thisDefaults{idxDefined(1)};
       if isnumeric(thisVal) && numel(thisVal) == 1
           thisDefaults(idxUndefined) = {NaN};
       elseif isnumeric(thisVal)
           thisDefaults(idxUndefined) = {'not applicable'};
       elseif isa(thisVal, 'function_handle')
           thisDefaults(idxUndefined) = {'not applicable'};
       elseif ischar(thisVal)
           thisDefaults(idxUndefined) = {'not applicable'};
       elseif iscell(thisVal) && isnumeric(thisVal{1}) && numel(thisVal{1}) == 1
           thisDefaults(idxUndefined) = {{NaN}};
       elseif iscell(thisVal) 
           thisDefaults(idxUndefined) = {{'not applicable'}};
       else
           thisDefaults(idxUndefined) = {'not applicable'};
       end
    end
    
    % Set the field
    OS.NonConstantDefaults.(OS.NonConstantDefaultFields{i}) = thisDefaults;
end
