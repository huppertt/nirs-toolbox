function hF = createhdlfilter(this)
%

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.df1tsos;

this.sethdl_abstractsos(hF);

hF.NumStateSLtype = conv2sltype(this.filterquantizer, 'StateWordLength', 'NumStateFracLength', true);
hF.DenStateSLtype = conv2sltype(this.filterquantizer, 'StateWordLength', 'DenStateFracLength', true);


hF.SectionInputSLtype = conv2sltype(this.filterquantizer, 'SectionInputWordLength', 'SectionInputFracLength', true);
hF.SectionOutputSLtype = conv2sltype(this.filterquantizer, 'SectionOutputWordLength', 'SectionOutputfracLength', true);
hF.MultiplicandSLtype = conv2sltype(this.filterquantizer, 'MultiplicandWordLength', 'MultiplicandfracLength', true);
% [EOF]
