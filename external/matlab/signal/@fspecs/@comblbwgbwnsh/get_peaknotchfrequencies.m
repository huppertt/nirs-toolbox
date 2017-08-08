function freqs = get_peaknotchfrequencies(this, freqs)
%GET_PEAKNOTCHFREQUENCIES PreGet function for the 'PeakNotchFrequencies'
%property

%   Copyright 2008 The MathWorks, Inc.

N =this.NumPeaksOrNotches;

if ~rem(N,2)
    freqs = (0:N/2)*(2/N);
else
    freqs = (0:(N-1)/2)*(2/N);
end

if ~this.NormalizedFrequency
    freqs = freqs*this.Fs/2;
end

% [EOF]
