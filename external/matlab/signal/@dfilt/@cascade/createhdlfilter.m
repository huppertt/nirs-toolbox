function hF = createhdlfilter(this)
%CREATHDLFILTER <short description>
%   OUT = CREATHDLFILTER(ARGS) <long description>

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.dfiltcascade;
sethdl_cascade(this, hF);
% [EOF]
