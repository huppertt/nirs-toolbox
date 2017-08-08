function freqs = get_peakfrequencies(this, freqs)
%GET_PEAKFREQUENCIES PreGet function for the 'PeakFrequencies'
%property

%   Copyright 2008 The MathWorks, Inc.

if isequal(lower(this.CombType),'peak')
    freqs = get(this.CurrentSpecs, 'PeakNotchFrequencies');
else
    freqs = ('N/A (notching comb)');
end

% [EOF]
