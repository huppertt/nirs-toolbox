function flag = ishalfnyqinterval(this)
%ISHALFNYQINTERVAL   True if the frequency response was calculated for only
%                    half the Nyquist interval.

%   Author(s): P. Pacheco
%   Copyright 1988-2003 The MathWorks, Inc.

flag = false;
if strcmpi(get(this,getrangepropname(this)),'half'),
    flag = true;
end

% [EOF]
