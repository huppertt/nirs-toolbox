function varargout = groupdelay(this, varargin)
%GROUPDELAY Group delay of a discrete-time filter.
%   Gd = GROUPDELAY(Hb) returns group delay response object Gd.
%
%   For additional parameters, see DFILT.BASEFILTER/FREQRESP.
%
%   See also DFILT, SIGNAL/GRPDELAY.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

hopts = uddpvparse('dspopts.freqresp', varargin{:});

inputs = freqzinputs(hopts);

[Gd, w] = base_resp(this, 'computegrpdelay', inputs{:});

opts = {};
if ~hopts.NormalizedFrequency
    opts = {'Fs', hopts.Fs};
end

if strcmpi(hopts.FrequencySpecification, 'NFFT')
    opts = {opts{:}, 'SpectrumRange', hopts.SpectrumRange};
end

h = dspdata.groupdelay(Gd, w, opts{:});

if nargout,
    varargout = {h};
else
    plot(h);
end

% [EOF]
