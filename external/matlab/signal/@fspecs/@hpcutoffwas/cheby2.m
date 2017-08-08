function Hd = cheby2(this, varargin)
%CHEBY2 Chebyshev Type II digital filter design.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

N = this.FilterOrder;

nfreq = get(this, 'NormalizedFrequency');
normalizefreq(this, true);

Fc = this.F3dB;

normalizefreq(this, nfreq);

Ast = this.Astop;

% Compute analog frequency
Wc = 1/tan(pi*Fc/2);

% Find corresponding analog stopband-edge frequency
Wst = Wc*cosh(1/N*acosh(sqrt(10^(Ast/10)-1)));

% Convert analog stopband-edge frequency to digital
Fst = 2*atan(1/Wst)/pi;

% Convert to highpass with stopband-edge specifications
hs = fspecs.hpstop(N,Fst,Ast);

Hd = cheby2(hs,varargin{:});

% [EOF]
