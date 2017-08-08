function [Pxx, Pxxc, f] = psd(varargin)
%PSD Power Spectral Density estimate.
%   PSD has been replaced by PWELCH.
%   Type help PWELCH for details.
%
%   See also PERIODOGRAM, PBURG, PCOV, PMCOV, PYULEAR, PMTM, PMUSIC, PEIG.

%   Author(s): T. Krauss, 3-26-93
%   Copyright 1988-2012 The MathWorks, Inc.

%   NOTE 1: To express the result of PSD, Pxx, in units of
%           Power per Hertz multiply Pxx by 1/Fs [1].
%
%   NOTE 2: The Power Spectral Density of a continuous-time signal,
%           Pss (watts/Hz), is proportional to the Power Spectral 
%           Density of the sampled discrete-time signal, Pxx, by Ts
%           (sampling period). [2] 
%       
%               Pss(w/Ts) = Pxx(w)*Ts,    |w| < pi; where w = 2*pi*f*Ts

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice hall, 1997, pg, 15
%     [2] A.V. Oppenheim and R.W. Schafer, Discrete-Time Signal
%         Processing, Prentice-Hall, 1989, pg. 731
%     [3] A.V. Oppenheim and R.W. Schafer, Digital Signal
%         Processing, Prentice-Hall, 1975, pg. 556

narginchk(1,7)

warning(message('signal:psd:PSDisObsolete'));

x = varargin{1};
[msg,nfft,Fs,window,noverlap,p,dflag,msgobj]=psdchk(varargin(2:end),x);
if ~isempty(msg), error(msgobj); end

% compute PSD
window = window(:);
n = length(x);		    % Number of data points
nwind = length(window); % length of window
if n < nwind            % zero-pad x if it has length less than the window length
    x(nwind)=0;  n=nwind;
end
% Make sure x is a column vector; do this AFTER the zero-padding 
% in case x is a scalar.
x = x(:);		

k = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
                    					% (k = fix(n/nwind) for noverlap=0)
%   if 0
%       disp(sprintf('   x        = (length %g)',length(x)))
%       disp(sprintf('   y        = (length %g)',length(y)))
%       disp(sprintf('   nfft     = %g',nfft))
%       disp(sprintf('   Fs       = %g',Fs))
%       disp(sprintf('   window   = (length %g)',length(window)))
%       disp(sprintf('   noverlap = %g',noverlap))
%       if ~isempty(p)
%           disp(sprintf('   p        = %g',p))
%       else
%           disp('   p        = undefined')
%       end
%       disp(sprintf('   dflag    = ''%s''',dflag))
%       disp('   --------')
%       disp(sprintf('   k        = %g',k))
%   end

index = 1:nwind;
KMU = k*norm(window)^2;	% Normalizing scale factor ==> asymptotically unbiased
% KMU = k*sum(window)^2;% alt. Nrmlzng scale factor ==> peaks are about right

Spec = zeros(nfft,1); % Spec2 = zeros(nfft,1);
for i=1:k
    if strcmp(dflag,'none')
        xw = window.*(x(index));
    elseif strcmp(dflag,'linear')
        xw = window.*detrend(x(index));
    else
        xw = window.*detrend(x(index),'constant');
    end
    index = index + (nwind - noverlap);
    Xx = abs(fft(xw,nfft)).^2;
    Spec = Spec + Xx;
%     Spec2 = Spec2 + abs(Xx).^2;
end

% Select first half
if ~any(any(imag(x)~=0)),   % if x is not complex
    if rem(nfft,2),    % nfft odd
        select = (1:(nfft+1)/2)';
    else
        select = (1:nfft/2+1)';
    end
    Spec = Spec(select);
%     Spec2 = Spec2(select);
%    Spec = 4*Spec(select);     % double the signal content - essentially
% folding over the negative frequencies onto the positive and adding.
%    Spec2 = 16*Spec2(select);
else
    select = (1:nfft)';
end
freq_vector = (select - 1)*Fs/nfft;

% find confidence interval if needed
if (nargout == 3) || ((nargout == 0) && ~isempty(p)),
    if isempty(p),
        p = .95;    % default
    end
    % Confidence interval from Kay, p. 76, eqn 4.16:
    % (first column is lower edge of conf int., 2nd col is upper edge)
    confid = Spec*chi2conf(p,k)/KMU;

    if noverlap > 0
        warning(message('signal:psd:InaccurateResult'));
    end
end

Spec = Spec*(1/KMU);   % normalize

% set up output parameters
if (nargout == 3),
   Pxx = Spec;
   Pxxc = confid;
   f = freq_vector;
elseif (nargout == 2),
   Pxx = Spec;
   Pxxc = freq_vector;
elseif (nargout == 1),
   Pxx = Spec;
elseif (nargout == 0),
   if ~isempty(p),
       P = [Spec confid];
   else
       P = Spec;
   end
   newplot;
   plot(freq_vector,10*log10(abs(P))), grid on
   xlabel('Frequency'), ylabel('Power Spectrum Magnitude (dB)');
end
