function flag = ishalfnyqinterval(this)
%ISHALFNYQINTERVAL   True if the spectrum was calculated for only half the
%                    Nyquist interval.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

flag = false;
if strcmpi(get(this,getrangepropname(this)),'onesided'),
    flag = true;
end

% [EOF]
