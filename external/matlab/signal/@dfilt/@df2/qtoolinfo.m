function info = qtoolinfo(this)
%QTOOLINFO   Return the information for the qtool.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

info = iir_qtoolinfo(this);

info.state.setops  = {'FracLabels', {'State'}};
info.state.syncops = [];

% [EOF]
