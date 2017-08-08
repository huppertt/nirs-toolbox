function Hd = cheby1(this, varargin)
%CHEBY1 Chebyshev Type I digital filter design.

%   Author(s): J. Schickler
%   Copyright 1988-2005 The MathWorks, Inc.

N = this.FilterOrder;

nfreq = get(this, 'NormalizedFrequency');
normalizefreq(this, true);

Fc1 = this.F3dB1;
Fc2 = this.F3dB2;
BWp = this.BWpass;

normalizefreq(this, nfreq);

[isvalid, errmsg, errid] = checkincfreqs(this,{'F3dB1','F3dB2'});
if ~isvalid
    error(message(errid,errmsg));
end

% If design comes from a filter designer object, convert propNames to valid
% designer object names.
propNames = {'BWp','F3dB2','F3dB1'};
if isprop(this,'FromFilterDesigner') && this.FromFilterDesigner
  propNames = signal.internal.filterdesigner.convertpropnames(propNames);
elseif isprop(this,'FromDesignfilt') && this.FromDesignfilt
  propNames = signal.internal.DesignfiltProcessCheck.convertpropnames(propNames,'fdesignToDesignfilt');
end

if BWp > (Fc2-Fc1),
    error(message('signal:fspecs:bpcutoffwbwp:cheby1:invalidBWSpec', propNames{1}, [propNames{2} '-' propNames{3}]));
end

% Make up a 3dB point for digital lowpass (arbitrarily)
Fc = (Fc1+Fc2)/2; 

% Compute alpha
alpha = cos((Fc2+Fc1)*pi/2)/cos((Fc2-Fc1)*pi/2); 
% Compute k
k = tan(Fc*pi/2)/tan((Fc2-Fc1)*pi/2);

% Compute theoretical lowpass passband edge
Wc = tan(pi*Fc/2);

% Determine lowpass passband edge
Fp = 2*atan(k*tan((BWp)*pi/2))/pi;

% Determine analog lowpass passband edge
Wp = tan(pi*Fp/2);

% Find epass, Apass
ep = 1/cosh((N/2)*acosh(Wc/Wp));
Ap = 10*log10(ep^2+1);

c1 = 2*alpha*k/(k+1);
c2 = (k-1)/(k+1);

% Solve LP to BP inverse Constantinides
j = complex(0,1);
z = exp(j*Fp*pi);
b = c1*(z+1);
ax2 = 0.5/(c2*z+1);
b2_4ac = sqrt(c1^2*(1+z)^2-4*(c2*z+1)*(z+c2));
Fp1 = abs(real(log(ax2*(b + b2_4ac))/(j*pi)));
Fp2 = abs(real(log(ax2*(b - b2_4ac))/(j*pi)));

% Convert to bandpass with passband-edge specifications
hs = fspecs.bppass(N,Fp1,Fp2,Ap);

Hd = cheby1(hs,varargin{:});

% [EOF]
