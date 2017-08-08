function b = isloggingon(this)
%ISLOGGINGON   True if the filter logging is on.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.

f = fipref;
b =  strcmpi(f.LoggingMode, 'on') && ...
        any(strmatch(f.DataTypeOverride , {'ForceOff','ScaledDoubles'}));


% [EOF]
