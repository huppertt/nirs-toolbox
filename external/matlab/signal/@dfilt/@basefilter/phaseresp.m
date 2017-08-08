function varargout = phase(this, varargin)
%PHASE   Phase response of a discrete-time filter.
%   [Phi,W] = PHASE(Hb,N) returns length N vectors Phi and W containing
%   the phase response of the discrete-time filter Hb, and the frequencies
%   (in radians) at which it is evaluated. The phase response is evaluated
%   at N points equally spaced around the upper half of the unit circle. If
%   you don't specify N, it defaults to 8192.
%
%   [Phi,W] = PHASE(Hb) returns a matrix Phi if Hb is a vector.  Each
%   column of the matrix corresponds to each filter in the vector.  If a
%   row vector of frequency points is specified, each row of the matrix
%   corresponds to each filter in the vector.
%
%   For additional parameters, see SIGNAL/PHASEZ.
%
%   See also DFILT, SIGNAL/PHASEZ.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

hopts = uddpvparse('dspopts.freqresp', varargin{:});

inputs = freqzinputs(hopts);

% Strings are faster than Function Handles
[Ph, w] = base_resp(this, 'computephasez', inputs{:});

opts = {};
if ~hopts.NormalizedFrequency
    opts = {'Fs', hopts.Fs};
end

if strcmpi(hopts.FrequencySpecification, 'NFFT')
    opts = {opts{:}, 'SpectrumRange', hopts.SpectrumRange};
end

h = dspdata.phaseresp(Ph, w, opts{:});

if nargout,
    varargout = {h};
else
    plot(h);
end


% [EOF]
