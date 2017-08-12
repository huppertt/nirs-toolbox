function output = isoptimargdbl(caller, argnames, varargin)
%

%ISOPTIMARGDBL returns true when input arguments are double or single.
%
% This is a helper function.

%   Copyright 2012 The MathWorks, Inc.


% Assume all the inputs are double
TF = true;
% Now, perform all the nagative tests and break asap.
try
    % SUPERIORFLOAT errors when superior input is neither single nor double;
    % We use try-catch to override SUPERIORFLOAT's error message when input
    % data type is integer.
    if ~strcmp(superiorfloat(varargin{:}),'double')
        TF = false;
    end
catch ME %#ok superiorfloat errors when superior input is neither single nor double.
    TF = false;
end

% It is not sufficient to check just the superior data type. We also want
% to make sure that input arguments are numeric.
if TF
    % Assumption: inputs at the end are more likely to be error prone.
    for i = length(varargin):-1:1 
        if ~isempty(varargin{i}) && ~isa(varargin{i},'double')
            TF = false;
            break;
        end
    end
end

if ~TF
    % We know at least some of the data type is wrong. Build an error
    % string containing the name of bad inputs. 
    nondblArgs = '';
    for i = 1:length(varargin)
        if ~isempty(varargin{i}) && ~isa(varargin{i},'double')
            nondblArgs = [nondblArgs, '''', argnames{i}, '''', ',']; %#ok
        end
    end
    if ~isempty(nondblArgs) % Remove the last comma.
        nondblArgs(end) = '';
    end
    output = getString(message('optimlib:commonMsgs:NonDoubleInput',caller,nondblArgs));
else
    output = '';
end

