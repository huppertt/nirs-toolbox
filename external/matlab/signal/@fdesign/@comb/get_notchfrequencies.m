function freqs = get_notchfrequencies(this, freqs)
%GET_NOTCHFREQUENCIES PreGet function for the 'NotchFrequencies'
%property

%   Copyright 2008 The MathWorks, Inc.

if isequal(lower(this.CombType),'notch')
    freqs = get(this.CurrentSpecs, 'PeakNotchFrequencies');
else
    freqs = ('N/A (peaking comb)');
end

% [EOF]
