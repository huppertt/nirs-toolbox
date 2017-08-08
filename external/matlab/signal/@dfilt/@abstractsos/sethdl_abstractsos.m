function sethdl_abstractsos(this, hhdl)
%SETHDLPROPSBASEFILTER Set the common props for HDLFILTER  from filter
%object
%   OUT = SETHDL_ABSTRACTSOS(ARGS) <long description>

%   Copyright 2007 The MathWorks, Inc.

this.sethdl_abstractfilter(hhdl);
[hhdl.RoundMode, hhdl.OverflowMode] = conv2hdlroundoverflow(this);

coeffs = coefficients(this);
hhdl.Coefficients = coeffs{1};
hhdl.scaleValues = this.ScaleValues;

hhdl.SectionOrder = secorder(this);
hhdl.NumSections = nsections(this);

hhdl.ScaleSLType = conv2sltype(this.filterquantizer, 'CoeffWordlength', 'ScaleValueFracLength');

hhdl.NumCoeffSLType = conv2sltype(this.filterquantizer, 'CoeffWordlength', 'NumFraclength');
hhdl.DenCoeffSLType = conv2sltype(this.filterquantizer, 'CoeffWordlength', 'DenFraclength');

hhdl.NumProdSLType = conv2sltype(this.filterquantizer, 'ProductWordlength', 'NumProdFraclength', true);
hhdl.DenProdSLType = conv2sltype(this.filterquantizer, 'ProductWordlength', 'DenProdFraclength', true);

hhdl.NumAccumSLType = conv2sltype(this.filterquantizer, 'AccumWordlength', 'NumAccumFraclength', true);
hhdl.DenAccumSLType = conv2sltype(this.filterquantizer, 'AccumWordlength', 'DenAccumFraclength', true);

% [EOF]
