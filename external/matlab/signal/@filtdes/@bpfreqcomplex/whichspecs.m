function specs = whichspecs(h)
%WHICHSPECS Determine which specs are required for this class.

%   Author(s): J. Schickler
%   Copyright 1988-2002 The MathWorks, Inc.

% Prop name, data type, default value, listener callback
specs(1) = cell2struct({'Fstop1','double',-16800,[],'freqspec'},specfields(h),2);
specs(2) = cell2struct({'Fpass1','double',-14400,[],'freqspec'},specfields(h),2);
specs(3) = cell2struct({'Fpass2','double',-12000,[],'freqspec'},specfields(h),2);
specs(4) = cell2struct({'Fstop2','double',-9600,[],'freqspec'},specfields(h),2);
specs(5) = cell2struct({'Fstop3','udouble',7200,[],'freqspec'},specfields(h),2);
specs(6) = cell2struct({'Fpass3','udouble',9600,[],'freqspec'},specfields(h),2);
specs(7) = cell2struct({'Fpass4','udouble',12000,[],'freqspec'},specfields(h),2);
specs(8) = cell2struct({'Fstop4','udouble',14400,[],'freqspec'},specfields(h),2);

% [EOF]
