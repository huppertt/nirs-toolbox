function hF = createhdlfilter(this)
%CREATEHDLFILTER Returns the corresponding hdlfiltercomp for HDL Code
%generation.

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.dffirt;

this.sethdl_dtffir(hF);

hF.StateSLtype = conv2sltype(this.filterquantizer, 'StateWordLength', 'StateFracLength', true);

% [EOF]
