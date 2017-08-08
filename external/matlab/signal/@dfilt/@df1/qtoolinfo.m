function info = qtoolinfo(this)
%QTOOLINFO   Return the information needed by the qtool.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

info = iir_qtoolinfo(this);

info.output.setops  = {'AutoscaleAvailable', 'Off'};
info.output.syncops = [];

% [EOF]
