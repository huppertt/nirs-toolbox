function privfq = set_privfq(this, privfq)
%SET_PRIVFQ   PreSet function for the 'privfq' property.

%   Author(s): J. Schickler
%   Copyright 1988-2006 The MathWorks, Inc.

% Recreate the listeners on the private filter quantizers.
l  = [handle.listener(privfq, 'quantizecoeffs', @super_quantizecoeffs); ...
    handle.listener(privfq, 'quantizestates', @lcl_quantizestates);...
    handle.listener(privfq, 'quantizefracdelay', @lcl_quantizefracdelay)];

set(l,  'callbacktarget', this);

set(this, 'filterquantizerlisteners', l);

% -------------------------------------------------------------------------
function lcl_quantizestates(this, eventData)

quantizestates(this);

% -------------------------------------------------------------------------
function lcl_quantizefracdelay(this, eventData)

quantizefd(this);

% [EOF]
