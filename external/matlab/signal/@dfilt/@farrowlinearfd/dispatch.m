function Hd = dispatch(this)
%DISPATCH   

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

d = this.FracDelay;
refd = this.reffracdelay;
for i=1:length(d),
    Hd(i) = lwdfilt.tf([1-d(i) d(i)]);
    Hd(i).refNum = [1-refd(i) refd(i)];
end


% [EOF]
