function s = whichspecs(h)
%WHICHSPECS

%   Author(s): J. Schickler
%   Copyright 1984-2003 The MathWorks, Inc.

s = [cell2struct({'Fnotch','udouble',9600,[],'freqspec'},specfields(h),2) fin_whichspecs(h)];

% [EOF]
