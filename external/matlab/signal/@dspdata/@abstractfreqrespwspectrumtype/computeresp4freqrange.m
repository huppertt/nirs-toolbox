function [H,W] = computeresp4freqrange(this,ishalfrange,ispsd,isnormalized,isDCcentered)
%COMPUTERESP4FREQRANGE   Compute the spectrum over the frequency range
%                        requested by the user.

%   Author(s): P. Pacheco
%   Copyright 1988-2004 The MathWorks, Inc.

% Define a boolean flag representing the state of SpectrumType property.
isonesided = ishalfnyqinterval(this);

% Make sure that Fs, frequency, and NormalizedFrequency property are all
% consistent.
normalizefreq(this,logical(isnormalized));

if ~ishalfrange && isonesided,    % User requested 'twosided' but obj has 'onesided'
    twosided(this);

elseif ishalfrange && ~isonesided, % User requested 'onesided' but obj has 'twosided'
    onesided(this);
end

if isDCcentered,
    centerdc(this);
end

H = this.Data;
W = this.Frequencies;

% Allow concrete classes to do further calculations if necessary.
[H,W] = thiscomputeresp4freqrange(this,H,W,ispsd);

% [EOF]
