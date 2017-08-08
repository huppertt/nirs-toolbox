function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Call alternate method
specs = mf_whichspecs(h);
specs(end+1) = cell2struct({'Apass','udouble',3,[],'magspec'},specfields(h),2);

% [EOF]
