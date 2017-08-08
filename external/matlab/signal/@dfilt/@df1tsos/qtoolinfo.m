function info = qtoolinfo(this)
%QTOOLINFO   Return the information needed by the QTool.

%   Author(s): J. Schickler
%   Copyright 1988-2009 The MathWorks, Inc.

info.coeff.setops  = {'Name', 'Coefficient', ...
    'FracLabels', {'Numerator', 'Denominator', 'Scale Values'}};
info.coeff.syncops = {'Num', 'Den', 'ScaleValue'};

info.product.setops  = {'FracLabels', {'Num.', 'Den.'}};
info.product.syncops = {'NumProd', 'DenProd'};

info.accum.setops  = {'FracLabels', {'Num.', 'Den.'}};
info.accum.syncops = {'NumAccum', 'DenAccum'};

info.state.setops  = {'FracLabels', {'Num.', 'Den.'}, 'AutoScaleAvailable', 'on'};
info.state.syncops = {'NumState', 'DenState'};

info.multiplicand.setops  = {};
info.multiplicand.syncops = {};

info.sectioninput.setops  = {'AutoScaleAvailable','on'};
info.sectioninput.syncops = {};

info.sectionoutput.setops  = {'AutoScaleAvailable','on'};
info.sectionoutput.syncops = {};

% [EOF]
