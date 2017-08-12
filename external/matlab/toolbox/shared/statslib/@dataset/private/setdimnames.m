function a = setdimnames(a,newnames)
%SETDIMNAMES Set dataset array DimNames property.

%   Copyright 2006-2011 The MathWorks, Inc.


if nargin < 2
    error(message('stats:dataset:setdimnames:TooFewInputs'));
end

if nargin == 2
    if numel(newnames) ~= a.ndims
        error(message('stats:dataset:setdimnames:WrongLength'));
    elseif ~iscell(newnames)
        error(message('stats:dataset:setdimnames:InvalidDimnames'));
    elseif ~all(cellfun(@isNonEmptyString,newnames))
        error(message('stats:dataset:setdimnames:InvalidDimnames'));
    end
    checkduplicatenames(newnames,'dimnames');
    a.props.DimNames = newnames(:)'; % a row vector
end

function tf = isNonEmptyString(s) % require a nonempty row of chars
tf = ischar(s) && isrow(s) && any(s ~= ' ');
