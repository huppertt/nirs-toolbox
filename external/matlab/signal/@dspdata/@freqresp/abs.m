function h = abs(this)
%ABS   Convert the frequency response to a magnitude response.

%   Author(s): J. Schickler
%   Copyright 2004 The MathWorks, Inc.

opts = {'SpectrumRange', this.SpectrumRange};
if ~this.NormalizedFrequency
    opts = {'Fs', this.Fs, opts{:}};
end

h = dspdata.magresp(this.Frequencies, abs(this.Data), opts{:});

% [EOF]
