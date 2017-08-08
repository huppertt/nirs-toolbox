function [Txy, f] = tfe(varargin)
%TFE Transfer Function Estimate.
%   TFE has been replaced by TFESTIMATE.  TFE still works but may be
%   removed in the future. Use TFESTIMATE instead.
%
%   TFESTIMATE does not support detrending.  Please use the DETREND
%   function if you need to detrend your signal. Type "help detrend" for
%   details.
%
%   See also PWELCH, CPSD, MSCOHERE,
%   ETFE, SPA, and ARX in the Identification Toolbox.

% 	Author(s): T. Krauss, 3-31-93
%   Copyright 1988-2004 The MathWorks, Inc.

error(nargchk(2,7,nargin,'struct'))

warning(message('signal:tfe:TFEisObsolete'));

x = varargin{1};
y = varargin{2};
[msg,nfft,Fs,window,noverlap,p,dflag,msgobj]=psdchk(varargin(3:end),x,y); %#ok
if ~isempty(msg), error(msgobj); end
    
% compute PSD and CSD
window = window(:);
n = length(x);		% Number of data points
nwind = length(window); % length of window
if n < nwind    % zero-pad x , y if length is less than the window length
    x(nwind)=0;
    y(nwind)=0;  
    n=nwind;
end
x = x(:);		% Make sure x is a column vector
y = y(:);		% Make sure y is a column vector
k = fix((n-noverlap)/(nwind-noverlap));	% Number of windows
					% (k = fix(n/nwind) for noverlap=0)
index = 1:nwind;

Pxx = zeros(nfft,1); Pxx2 = zeros(nfft,1);
Pxy = zeros(nfft,1); Pxy2 = zeros(nfft,1);
for i=1:k
    if strcmp(dflag,'none')
        xw = window.*x(index);
        yw = window.*y(index);
    elseif strcmp(dflag,'linear')
        xw = window.*detrend(x(index));
        yw = window.*detrend(y(index));
    else
        xw = window.*detrend(x(index),0);
        yw = window.*detrend(y(index),0);
    end
    index = index + (nwind - noverlap);
    Xx = fft(xw,nfft);
    Yy = fft(yw,nfft);
    Xx2 = abs(Xx).^2;
    Xy2 = Yy.*conj(Xx);
    Pxx = Pxx + Xx2;
    Pxx2 = Pxx2 + abs(Xx2).^2;
    Pxy = Pxy + Xy2;
    Pxy2 = Pxy2 + Xy2.*conj(Xy2);
end

% Select first half
if ~any(any(imag([x y])~=0)),   % if x and y are not complex
    if rem(nfft,2),    % nfft odd
        select = 1:(nfft+1)/2;
    else
        select = 1:nfft/2+1;   % include DC AND Nyquist
    end
    Pxx = Pxx(select);
    Pxy = Pxy(select);
else
    select = 1:nfft;
end
Trans = Pxy ./ Pxx;             % transfer function estimate 
freq_vector = (select - 1)'*Fs/nfft;

% set up output parameters
if (nargout == 2),
   Txy = Trans;
   f = freq_vector;
elseif (nargout == 1),
   Txy = Trans;
elseif (nargout == 0),   % do a plot
   newplot;
   plot(freq_vector,20*log10(abs(Trans))), grid on
   xlabel('Frequency'), ylabel('Transfer Function Estimate (dB)');
end
