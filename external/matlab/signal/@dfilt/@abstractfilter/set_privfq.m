function privfq = set_privfq(this, privfq)
%SET_PRIVFQ   PreSet function for the 'privfq' property.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

% Recreate the listeners on the private filter quantizers.
l  = [handle.listener(privfq, 'quantizecoeffs', @super_quantizecoeffs); ...
    handle.listener(privfq, 'quantizestates', @lcl_quantizestates)];

set(l,  'callbacktarget', this);

set(this, 'filterquantizerlisteners', l);

% -------------------------------------------------------------------------
function lcl_quantizestates(this, eventData)

quantizestates(this);

% [EOF]
