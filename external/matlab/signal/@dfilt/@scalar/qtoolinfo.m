function info = qtoolinfo(this)
%QTOOLINFO   Return the information for the qtool.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

info.coeff.setops = {'Name', 'Coefficient', 'FracLabels', {'Coefficient'}};
info.coeff.syncops = {'Coeff'};

info.accum   = [];
info.product = [];

% [EOF]
