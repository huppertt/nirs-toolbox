function Hd = cheby1(this, varargin)
%CHEBY1 Chebyshev Type I digital filter design.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

N = this.FilterOrder;

nfreq = get(this, 'NormalizedFrequency');
normalizefreq(this, true);

Fc = this.F3dB;

normalizefreq(this, nfreq);

Ap = this.Apass;

% Compute analog frequency
Wc = tan(pi*Fc/2);

% Find corresponding analog passband-edge frequency
Wp = Wc/cosh(1/N*acosh(1/sqrt(10^(Ap/10)-1)));

% Convert analog passband-edge frequency to digital
Fp = 2*atan(Wp)/pi;

% Convert to lowpass with passband-edge specifications
hs = fspecs.lppass(N,Fp,Ap);

Hd = cheby1(hs,varargin{:});

% [EOF]
