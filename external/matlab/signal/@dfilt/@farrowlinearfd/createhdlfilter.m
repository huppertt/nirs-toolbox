function hF = createhdlfilter(this)
%CREATEHDLFILTER Returns the corresponding hdlfiltercomp for HDL Code
%generation.

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.farrowlinearfd;
this.sethdl_abstractfarrow(hF);

hF.TapSumSLType = conv2sltype(this.filterquantizer, 'TapSumWordlength', 'TapSumFraclength', true);

% [EOF]
