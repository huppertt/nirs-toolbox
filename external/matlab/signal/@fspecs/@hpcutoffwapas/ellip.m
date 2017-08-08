function Hd = ellip(this, varargin)
%ELLIP Elliptic digital filter design.

%   Author(s): V. Pellissier
%   Copyright 2005-2006 The MathWorks, Inc.

nfreq = get(this, 'NormalizedFrequency');
normalizefreq(this, true);

% Design prototype filter for any fp, to keep things simple, we make fp
% equal to F3dB
hs = fspecs.lppassastop(this.FilterOrder,this.F3dB,this.Apass,this.Astop);

Hproto = ellip(hs,varargin{:});
Hproto.setfdesign(hs);

% Find 3 db point
measurements = measure(Hproto);

% Apply frequency transformation to obtain final filter 
Hd = iirlp2hp(Hproto,measurements.F3dB,this.F3dB);

normalizefreq(this, nfreq);

Hd.setfmethod(Hproto.getfmethod);

sosscale(Hproto.getfmethod, Hd);

% [EOF]