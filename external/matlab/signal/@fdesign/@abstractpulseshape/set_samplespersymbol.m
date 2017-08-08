function samplesPerSymbol = set_samplespersymbol(this, samplesPerSymbol)
%SET_SAMPLESPERSYMBOL PreSet function for the 'SamplesPerSymbol' property
%   OUT = SET_SAMPLESPERSYMBOL(ARGS) <long description>

%   Copyright 2008 The MathWorks, Inc.

set(this, 'privSamplesPerSymbol', samplesPerSymbol);

% If the SamplesPerSymbol property exists on the new specifications, set it.
if isprop(this.CurrentSpecs, 'SamplesPerSymbol')
    set(this.CurrentSpecs, 'SamplesPerSymbol', samplesPerSymbol);
end

% [EOF]
