function nfft = set_nfft(this, nfft) %#ok
%SET_NFFT   PreSet function for the 'nfft' property.

%   Author(s): P. Pacheco
%   Copyright 1988-2011 The MathWorks, Inc.

% Welch uses segment length instead of input length.
% nextpow2 = max(256,nextpow2(inputlength))
% auto = max(256,inputlength)

validStrs = {'Auto','Nextpow2'};

if isnumeric(nfft),
    if ~isscalar(nfft) || nfft<=0 || rem(nfft,1),
        error(message('signal:dspopts:abstractspectrumwfreqpoints:set_nfft:invalidNFFTValue', 'NFFT', 'NFFT', validStrs{ 1 }, validStrs{ 2 }));
    end

else
    idx = [];
    for k=1:length(validStrs),
        if regexp(lower(validStrs{k}),['^',lower(nfft)],'once');
            idx=k;
        end
    end
    if isempty(idx),
        error(message('signal:dspopts:abstractspectrumwfreqpoints:set_nfft:invalidNFFTValue', 'NFFT', 'NFFT', validStrs{ 1 }, validStrs{ 2 }));
    else
        % Use full string with correct capitalization.
        nfft = validStrs{idx};
    end
end

% [EOF]
