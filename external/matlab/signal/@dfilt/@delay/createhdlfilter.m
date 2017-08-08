function hF = createhdlfilter(this)
%CREATHDLFILTER <short description>
%   OUT = CREATHDLFILTER(ARGS) <long description>

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.delay;
this.sethdl_abstractfilter(hF);
hF.Latency = this.Latency;
% [EOF]
