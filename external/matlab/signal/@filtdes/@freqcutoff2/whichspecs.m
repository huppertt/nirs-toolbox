function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): R. Losada
%   Copyright 1988-2002 The MathWorks, Inc.

% Prop name, data type, default value, listener callback
specs(1) = cell2struct({'Fc1','udouble',8400,[],'freqspec'},specfields(h),2);

specs(2) = cell2struct({'Fc2','udouble',13200,[],'freqspec'},specfields(h),2);


