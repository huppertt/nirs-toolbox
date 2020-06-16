function t = setDescription(t,newdescr)
%SETDESCRIPTION Set table Description property.

%   Copyright 2012 The MathWorks, Inc.

import matlab.internal.tableUtils.isstring

if ~isstring(newdescr)
    error(message('MATLAB:table:InvalidDescription'));
elseif isempty(newdescr)
    t.props.Description = '';
else
    t.props.Description = newdescr;
end
