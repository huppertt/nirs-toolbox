function s = whichspecs(h)
%WHICHSPECS

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

s = [cell2struct({'Fpeak','udouble',9600,[],'freqspec'},specfields(h),2) fin_whichspecs(h)];

% [EOF]
