function type = set_type(this, type)
%SET_TYPE PreSet function for the 'type' property

%   Copyright 2007 The MathWorks, Inc.

this.privType = type;
updatecurrentspecs(this)

% [EOF]
