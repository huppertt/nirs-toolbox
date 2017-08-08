function info = qtoolinfo(this)
%QTOOLINFO   Return the information for the qtool.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

info.coeff.setops = {'Name', 'Lattice', ...
    'FracLabels', {'Lattice'}};
info.coeff.syncops = {'Lattice'};

info.state.setops  = {'FracLabels', {'State'}};
info.state.syncops = [];

info.product.setops  = {'FracLabels', {'Product'}};
info.product.syncops = {'Product'}; % Syncs the defaults.

info.accum.setops  = {'FracLabels', {'Accum.'}};
info.accum.syncops = {'Accum'};

% [EOF]
