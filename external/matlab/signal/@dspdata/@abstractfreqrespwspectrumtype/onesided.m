function onesided(this)
%ONESIDED   Convert a two-sided spectrum to a one-sided spectrum.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

newSpectrumType = 'onesided';
if strcmpi(this.SpectrumType,newSpectrumType),
    return;     % Spectrum already one-sided.
end

% Force data and frequencies to be in the range 0-Fs.
centerdc(this,false);  % no-op if dc is not centered.

Pxx = this.Data;
W   = this.Frequencies;
[nfft,nchans] = size(Pxx);

% Convert a 'twosided' spectrum (and frequencies) to a 'onesided' spectrum.
if rem(nfft,2),
    select = 1:(nfft+1)/2;                    % ODD;  take only [0,pi)
    Pxx = [Pxx(1,:); 2*Pxx(select(2:end),:)]; % Don't multiply DC term by 2.
else
    select = 1:nfft/2+1;                      % EVEN; take only [0,pi]
    Pxx = [Pxx(1,:); 2*Pxx(select(2:end-1),:); Pxx(select(end),:)]; % Don't multiple DC & Nyquist by 2.
end
W = W(select);

this.Data = Pxx;
this.Frequencies = W;
setspectrumtype(this,newSpectrumType);

% [EOF]
