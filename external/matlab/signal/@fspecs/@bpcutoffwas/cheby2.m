function Hd = cheby2(this, varargin)
%CHEBY2 Chebyshev Type II digital filter design.

%   Copyright 1988-2012 The MathWorks, Inc.

N   = this.FilterOrder;
Ast = this.Astop;

nfreq = get(this, 'NormalizedFrequency');
normalizefreq(this, true);

Fc1 = this.F3dB1;
Fc2 = this.F3dB2;

normalizefreq(this, nfreq);

[isvalid, errmsg, errid] = checkincfreqs(this,{'F3dB1','F3dB2'});
if ~isvalid
    error(message(errid,errmsg));
end

% Make up a 3dB point for digital lowpass (arbitrarily)
Fc = (Fc1+Fc2)/2; 

% Compute theoretical lowpass passband edge
Wc = tan(pi*Fc/2);

% Find corresponding analog stopband-edge frequency
Wst = Wc*cosh(1/(N/2)*acosh(sqrt(10^(Ast/10)-1))); % Use half the order

% Convert analog stopband-edge frequency to digital
Fst = 2*atan(Wst)/pi;

% Compute alpha
alpha = cos((Fc2+Fc1)*pi/2)/cos((Fc2-Fc1)*pi/2); 
% Compute k
k = tan(Fc*pi/2)/tan((Fc2-Fc1)*pi/2);

c1 = 2*alpha*k/(k+1);
c2 = (k-1)/(k+1);

% Solve LP to BP inverse Constantinides (Antoniou p.243)
j = complex(0,1);
z = exp(j*Fst*pi);
b = c1*(z+1);
ax2 = 0.5/(c2*z+1);
b2_4ac = sqrt(c1^2*(1+z)^2-4*(c2*z+1)*(z+c2));
Fst1 = abs(real(log(ax2*(b + b2_4ac))/(j*pi)));
Fst2 = abs(real(log(ax2*(b - b2_4ac))/(j*pi)));

% Convert to bandpass with stopband-edge specifications
hs = fspecs.bpstop(N,Fst1,Fst2,Ast);

Hd = cheby2(hs,varargin{:});

