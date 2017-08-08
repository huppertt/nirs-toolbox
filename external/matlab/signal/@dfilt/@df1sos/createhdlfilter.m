function hF = createhdlfilter(this)
%

%   Copyright 2007 The MathWorks, Inc.

hF = hdlfilter.df1sos;

this.sethdl_abstractsos(hF);

hF.NumStateSLtype = conv2sltype(this.filterquantizer, 'NumStateWordLength', 'NumStateFracLength');
hF.DenStateSLtype = conv2sltype(this.filterquantizer, 'DenStateWordLength', 'DenStateFracLength');

% [EOF]
