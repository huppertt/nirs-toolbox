function halfrange(this)
%HALFRANGE   Power spectrum calculated over half the Nyquist interval.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

newSpectrumRange = 'half';
if strcmpi(this.SpectrumRange,newSpectrumRange),
    % Spectrum already 'half' the Nyquist interval.
    return
end

% Convert a spectrum calculated over the 'whole' Nyquist interval to
% spectrum calculated over 'half' the Nquist interval.
[nfft,nchans] = size(this.Data);
if rem(nfft,2),    select = 1:(nfft+1)/2;  % ODD;  take only [0,pi)
else               select = 1:nfft/2+1;    % EVEN; take only [0,pi]
end

% Update object with new values.
set(this,...
    'Data',this.Data(select,:),...
    'Frequencies', this.Frequencies(select));

setspectrumtype(this,newSpectrumRange);

% [EOF]
