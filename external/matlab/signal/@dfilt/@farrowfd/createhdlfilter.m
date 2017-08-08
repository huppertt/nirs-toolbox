function hF = createhdlfilter(this)
%CREATEHDLFILTER Returns the corresponding hdlfiltercomp for HDL Code
%generation.

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.farrowfd;
this.sethdl_abstractfarrow(hF);

coeffs = coefficients(this);
hF.Coefficients = coeffs{1};

hF.CoeffSLType = conv2sltype(this.filterquantizer, 'CoeffWordlength', 'CoeffFraclength');
hF.MultiplicandSLType = conv2sltype(this.filterquantizer, 'MultiplicandWordlength', 'MultiplicandFraclength', true);
hF.FDProdSLType = conv2sltype(this.filterquantizer, 'FDProdWordlength', 'FDProdFraclength', true);

% [EOF]
