function freqresp = set_freqresp(this, freqresp)
%SET_FREQRESP   PreSet function for the 'freqresp' property.

%   Author(s): V. Pellissier
%   Copyright 2005 The MathWorks, Inc.


% Force row vector
freqresp = freqresp(:).';

% [EOF]
