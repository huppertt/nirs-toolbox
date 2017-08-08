function response = set_response(this, response)
%SET_RESPONSE   PreSet function for the 'Response' property.

%   Author(s): J. Schickler
%   Copyright 2005 The MathWorks, Inc.

set(this, 'privResponse', response);

% Update the currently stored FDESIGN object.
updatecurrentfdesign(this);

% [EOF]
