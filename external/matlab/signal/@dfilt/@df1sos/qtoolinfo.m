function info = qtoolinfo(this)
%QTOOLINFO   Return the information needed by the QTool.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

info.coeff.setops  = {'Name', 'Coefficient', ...
    'FracLabels', {'Numerator', 'Denominator', 'Scale Values'}};
info.coeff.syncops = {'Num', 'Den', 'ScaleValue'};

info.product.setops  = {'FracLabels', {'Num.', 'Den.'}};
info.product.syncops = {'NumProd', 'DenProd'};

info.accum.setops  = {'FracLabels', {'Num.', 'Den.'}};
info.accum.syncops = {'NumAccum', 'DenAccum'};

info.state.setops  = {'FracLabels', {'Num.', 'Den.'}, ...
    'Name', 'Num. state', 'WordLabel2', 'Den. state', 'AutoScaleAvailable', 'off'};
info.state.syncops = {'NumState', 'DenState'};

% There is no multiplicand, stage input, or stage output.

% [EOF]
