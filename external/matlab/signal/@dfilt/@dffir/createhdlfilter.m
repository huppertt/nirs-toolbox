function hF = createhdlfilter(this)
%CREATEHDLFILTER Returns the corresponding hdlfiltercomp for HDL Code
%generation.

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.dffir;
this.sethdl_dtffir(hF);
% [EOF]
