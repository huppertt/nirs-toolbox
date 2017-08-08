function varargout = phasedelay(this, varargin)
%PHASEDELAY   Compute the phase delay.
%   [PHI,W] = PHASEDELAY(B,A,N) returns the N-point phase delay response
%   vector PHI and the N-point frequency vector W in radians/sample of
%   the filter:
%               jw               -jw              -jmw 
%        jw  B(e)    b(1) + b(2)e + .... + b(m+1)e
%     H(e) = ---- = ------------------------------------
%               jw               -jw              -jnw
%            A(e)    a(1) + a(2)e + .... + a(n+1)e
%   given numerator and denominator coefficients in vectors B and A. The
%   phase response is evaluated at N points equally spaced around the
%   upper half of the unit circle. If N isn't specified, it defaults to
%   512.

%   Author(s): J. Schickler
%   Copyright 1988-2004 The MathWorks, Inc.

if nargin > 1 && ischar(varargin{1}) && ~any(strcmpi(varargin{1}, {'whole', 'half'}))
    newobj = true;
else
    newobj = false;
end

if newobj
    hopts = uddpvparse('dspopts.freqresp', varargin{:});
    
    inputs = hopts.freqzinputs;

    [Phi, W] = base_resp(this, 'computephasedelay', inputs{:});

    opts = {};
    if ~hopts.NormalizedFrequency
        opts = {'Fs', hopts.Fs};
    end

    if strcmpi(hopts.FrequencySpecification, 'NFFT')
        opts = {opts{:}, 'SpectrumRange', hopts.SpectrumRange};
    end

    h = dspdata.phasedelay(Phi, W, opts{:});
    
    if nargout
        varargout = {h};
    else
        plot(h);
    end
elseif nargout

    [Phi, W] = base_resp(this, 'computephasedelay', varargin{:});
    varargout = {Phi, W};
else
    [this, opts] = freqzparse(this, varargin{:});
    fvtool(this, 'phasedelay', opts);
end

% [EOF]
