function a = setdescription(a,newdescr)
%SETDESCRIPTION Set dataset array Description property.

%   Copyright 2006-2011 The MathWorks, Inc.


if nargin < 2
    error(message('stats:dataset:setdescription:TooFewInputs'));
end

if nargin == 2
    if isempty(newdescr)
        a.props.Description = '';
        return
    end
    a.props.Description = newdescr;
end
