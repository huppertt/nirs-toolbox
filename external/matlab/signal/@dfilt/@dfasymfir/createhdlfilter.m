function hF = createhdlfilter(this)
%CREATEHDLFILTER Returns the corresponding hdlfiltercomp for HDL Code
%generation.

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.dfasymfir;

this.sethdl_dtffir(hF);

hF.TapsumSLtype = conv2sltype(this.filterquantizer, 'TapsumWordLength', 'TapSumFracLength', true);

% [EOF]
