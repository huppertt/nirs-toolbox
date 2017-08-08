function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

% Prop name, data type, default value, listener callback
specs(1) = cell2struct({'DstopUpper','udouble',.01,[],'magspec'},specfields(h),2);
specs(2) = cell2struct({'DstopLower','udouble',.01,[],'magspec'},specfields(h),2);
specs(3) = cell2struct({'DpassUpper','udouble',.02,[],'magspec'},specfields(h),2);
specs(4) = cell2struct({'DpassLower','udouble',.02,[],'magspec'},specfields(h),2);

% [EOF]

