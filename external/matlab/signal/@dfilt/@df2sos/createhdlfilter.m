function hF = createhdlfilter(this)
%

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.df2sos;

this.sethdl_abstractsos(hF);

hF.SectionInputSLtype = conv2sltype(this.filterquantizer, 'SectionInputWordLength', 'SectionInputFracLength', true);
hF.SectionOutputSLtype = conv2sltype(this.filterquantizer, 'SectionOutputWordLength', 'SectionOutputfracLength', true);
hF.StateSLtype = conv2sltype(this.filterquantizer, 'StateWordLength', 'StateFracLength');

% [EOF]
