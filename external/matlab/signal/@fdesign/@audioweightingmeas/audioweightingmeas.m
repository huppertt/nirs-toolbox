function this = audioweightingmeas(hfilter, varargin)
%AUDIOWEIGHTINGMEAS   Construct a AUDIOWEIGHTINGMEAS object.

%   Copyright 2009 The MathWorks, Inc.

error(nargchk(1,inf,nargin,'struct'));

this = fdesign.audioweightingmeas;

% Parse the inputs.
minfo = parseconstructorinputs(this, hfilter, varargin{:});

% Measure attenuations at the frequencies specified by the standards
Resp = freqz(hfilter,minfo.F,minfo.Fs);
this.Magnitudes = 20*log10(abs(Resp));
this.Frequencies = minfo.F;

this.Magnitudes(this.Frequencies > minfo.Fs/2) = NaN;
this.Frequencies(this.Frequencies > minfo.Fs/2) = NaN;

% Measure attenuations using interpolated masks. These measurements will be the
% ones used to decide if specs are met or not.
Resp = freqz(hfilter,minfo.Finterp,minfo.Fs);
this.MagnitudesInterp = 20*log10(abs(Resp));
this.FrequenciesInterp = minfo.Finterp;

% [EOF]
