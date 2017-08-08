function disp(h)
%DISP   Object Display.

%   Author(s): P. Costa
%   Copyright 1988-2003 The MathWorks, Inc.

s = get(h);

if isfixed(h),
    [NumStr,DenStr] = getfistr(h);
    disp(changedisplay(s, 'Numerator', NumStr,'Denominator', DenStr))
else
    disp(s)
end


% -------------------------------------------------------------------------
function [NumStr,DenStr] = getfistr(h)

szNumS = size(h.Numerator);
szDenS = size(h.Denominator);
NumStr = ['[',num2str(szNumS(1)),'x',num2str(szNumS(2)),' ','fi',']'];
DenStr = ['[',num2str(szDenS(1)),'x',num2str(szDenS(2)),' ','fi',']'];

% [EOF]
