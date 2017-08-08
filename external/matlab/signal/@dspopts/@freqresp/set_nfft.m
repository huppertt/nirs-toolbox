function nfft = set_nfft(this, nfft)
%SET_NFFT   PreSet function for the 'nfft' property.

%   Author(s): J. Schickler
%   Copyright 2004-2011 The MathWorks, Inc.

set(this, 'FrequencySpecification', 'NFFT');

if ischar(nfft) || ~isscalar(nfft) || nfft<=0 || rem(nfft,1),
    error(message('signal:dspopts:freqresp:set_nfft:invalidDataType', 'NFFT', 'NFFT'));
end
% [EOF]
