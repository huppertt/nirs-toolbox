function info = qtoolinfo(this)
%QTOOLINFO   Return the information for the qtool.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

info = iir_qtoolinfo(this);

info.state.setops  = {'FracLabels', {'Num.', 'Den.'}, 'AutoScaleAvailable', 'On'};
info.state.syncops = {'NumState', 'DenState'};

info.multiplicand.setops  = [];
info.multiplicand.syncops = [];

% [EOF]
