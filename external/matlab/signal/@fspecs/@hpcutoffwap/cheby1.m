function Hd = cheby1(this, varargin)
%CHEBY1 Chebyshev Type I digital filter design.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

N = this.FilterOrder;
Ap = this.Apass;

nfreq = get(this, 'NormalizedFrequency');
normalizefreq(this, true);

Fc = this.F3db;

normalizefreq(this, nfreq);

% Compute analog frequency
Wc = 1/tan(pi*Fc/2);

% Find corresponding analog passband-edge frequency
Wp = Wc/cosh(1/N*acosh(1/sqrt(10^(Ap/10)-1)));

% Convert analog passband-edge frequency to digital
Fp = 2*atan(1/Wp)/pi;

% Convert to highpass with passband-edge specifications
hs = fspecs.hppass(N,Fp,Ap);

Hd = cheby1(hs,varargin{:});

% [EOF]
