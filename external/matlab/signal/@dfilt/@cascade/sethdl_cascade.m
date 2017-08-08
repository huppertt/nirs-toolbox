function sethdl_cascade(this, hhdl)
%SETHDL_CASCADE Set the properties of hdlfiltercomp (hhdl) from the
%filter object.

%   Copyright 2007 The MathWorks, Inc.

for n = 1: length(this.Stage)
    hFn = createhdlfilter(this.Stage(n));
    hhdl.Stage = [hhdl.Stage; hFn];
end
hhdl.RateChangeFactors = this.getratechangefactors;

% [EOF]
