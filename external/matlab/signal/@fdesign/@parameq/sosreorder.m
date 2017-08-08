function sosreorder(this,Hd)
%SOSREORDER   Reorder SOS filter.

%   Copyright 2008 The MathWorks, Inc.

if islphpreorder(this.CurrentSpecs),
    lphpreorder(this,Hd);
else
    bpbsreorder(this,Hd);
end

% [EOF]
