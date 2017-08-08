function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Call super's method
specs = mf_whichspecs(h);

% Prop name, data type, default value, listener callback
specs(end+1) = cell2struct({'Apass1','udouble',1,[],'magspec'},specfields(h),2);
specs(end+1) = cell2struct({'Wstop','udouble',1,[],''},specfields(h),2);
specs(end+1) = cell2struct({'Wpass2','udouble',1,[],''},specfields(h),2);

% [EOF]
