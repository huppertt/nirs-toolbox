function sethdl_dtffir(this, hhdl)
%SETHDLPROPSBASEFILTER Set the common props for HDLFILTER  from filter
%object

%   Copyright 2007 The MathWorks, Inc.
this.sethdl_abstractfilter(hhdl)
[hhdl.RoundMode, hhdl.OverflowMode] = conv2hdlroundoverflow(this);

coeffs = coefficients(this);
hhdl.Coefficients = coeffs{1};

hhdl.CoeffSLType = conv2sltype(this.filterquantizer, 'CoeffWordlength', 'NumFraclength');
hhdl.ProductSLType = conv2sltype(this.filterquantizer, 'ProductWordLength', 'ProductFracLength', true);
hhdl.AccumSLType = conv2sltype(this.filterquantizer, 'AccumWordLength', 'AccumFracLength', true);


% [EOF]
