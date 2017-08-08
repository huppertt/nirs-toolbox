function c = cparam(h)
%CPARAM   

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

F1 = h.Fstop1;
F2 = h.Fstop2;
c = computecparam(h,F1,F2);

% [EOF]
