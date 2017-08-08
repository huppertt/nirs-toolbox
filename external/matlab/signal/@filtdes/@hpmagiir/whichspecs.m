function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Call super's method
specs = mis_whichspecs(h);

% Prop name, data type, default value, listener callback
specs(end+1) = cell2struct({'Apass','udouble',1,[],'magspec'},specfields(h),2);






