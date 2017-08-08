function combtype = set_combtype(this, combtype)
%SET_COMBTYPE PreSet function for the 'CombType' property

%   Copyright 2008 The MathWorks, Inc.

%This function connects the set value of CombType in fdesign with the value
%in fspecs
set(this.CurrentSpecs, 'CombType', combtype);

% [EOF]
