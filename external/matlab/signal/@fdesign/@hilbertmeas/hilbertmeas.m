function this = hilbertmeas(hfilter, varargin)
%HILBERTMEAS   Construct a HILBERTMEAS object.

%   Author(s): P. Costa & J. Schickler
%   Copyright 2005 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

this = fdesign.hilbertmeas;

minfo = parseconstructorinputs(this, hfilter, varargin{:});

this.TransitionWidth = minfo.TransitionWidth;

% Apass represents the passband ripple!
if this.NormalizedFrequency, Fs = 2;
else,                        Fs = this.Fs; end

wpass1 = this.TransitionWidth/2;
wpass2 = Fs/2-this.TransitionWidth/2;

this.Apass = measureripple(this, hfilter, wpass1, wpass2, minfo.Apass);

% [EOF]
