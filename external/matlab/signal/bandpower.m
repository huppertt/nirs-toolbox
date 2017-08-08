function pwr = bandpower(varargin)
%BANDPOWER Band power
%   P = BANDPOWER(X) computes the average power in the input signal vector,
%   X.  If X is a matrix, then BANDPOWER computes the average power in each
%   column independently.
%
%   P = BANDPOWER(X, Fs, FREQRANGE) specifies FREQRANGE as a two-element
%   vector of real values, specifying the two frequencies between which you
%   want to measure the power.  Fs is the sampling rate of the input
%   signal.  The power is estimated by applying a Hamming window and using
%   a periodogram of the same length as the input vector.  If the input
%   vector X contains N samples, then FREQRANGE must be contained within
%   the interval:
%      [              0, Fs/2          ] if X is real and N is even
%      [              0, Fs*(N-1)/(2*N)] if X is real and N is odd
%      [-Fs*(N-2)/(2*N), Fs/2          ] if X is complex and N is even
%      [-Fs*(N-1)/(2*N), Fs*(N-1)/(2*N)] if X is complex and N is odd
%
%   P = BANDPOWER(Pxx, F, 'psd') computes the average power via a rectangle
%   approximation of the integral of the Power Spectral Density (PSD)
%   estimate, given in the vector Pxx over the frequencies specified in
%   vector F.
%
%   P = BANDPOWER(Pxx, F, FREQRANGE, 'psd') specifies FREQRANGE as a
%   two-element vector of real values, specifying the two frequencies
%   between which you want to measure the power.  F is a vector containing
%   the frequencies that correspond to the estimates given in Pxx.
%
%   NOTE: If the frequency range values don't exactly match the frequency
%   values stored in F then the next closest value is used.
%
%   % Example
%   %   Examine the power present within 50 kHz of the carrier of a 1 kHz
%   %   sinusoid FM modulated at 250 kHz with a modulation index of 0.1.
%
%   Fs = 1.0e6; Fc = 250e3; Fd = 1e3; m = 0.1;
%   x = vco(m*sin(2*pi*Fd*(1:1000)/Fs),Fc,Fs);
%   periodogram(x,[],length(x),Fs);
%   mod_power_dB = 10*log10(bandpower(x, Fs, Fc + [-50e3 50e3]))
%
%   See also POWERBW OBW PERIODOGRAM PWELCH MEANFREQ MEDFREQ PLOMB.

%   Copyright 2012-2014 The MathWorks, Inc.

narginchk(1,4);

matches = find(strcmpi('psd',varargin));
varargin(matches) = [];

if any(matches)
    pwr = psdbandpower(varargin{:});
else
    pwr = timedomainbandpower(varargin{:});
end

function pwr = timedomainbandpower(x, fs, freqrange)

% perform column vector conversion before checking 2d matrix
if isvector(x)
  x = x(:);
end

validateattributes(x, {'numeric'},{'2d','finite'}, ...
    'bandpower', 'x', 1);
  
if nargin==1
    % full range specified
    pwr = rms(x).^2;
    return
elseif nargin == 2
    error(message('signal:bandpower:FreqRangeMissing'));
else
    validateattributes(fs, {'numeric'},{'scalar','finite','real','positive'}, ...
        'bandpower', 'Fs', 2);
    % Compute periodogram using a hamming window with the same length as
    % input
    
    % Cast to enforce Precision rules
    fs = double(fs);
    freqrange = signal.internal.sigcasttofloat(freqrange,'double',...
      'bandpower','FREQRANGE','allownumeric');
    
    n = size(x,1);
    
    if isreal(x)
        [Pxx, F] = periodogram(x, hamming(n), n, fs);
    else
        [Pxx, F] = periodogram(x, hamming(n), n, fs, 'centered');
    end
    % Return the bandpower
    pwr = psdbandpower(Pxx, F, freqrange);
end

function pwr = psdbandpower(Pxx, W, freqrange)

% perform column vector conversion before checking 2d matrix
if isvector(Pxx)
  Pxx = Pxx(:);
end

validateattributes(Pxx, {'numeric'},{'2d','finite','real'}, ...
    'bandpower', 'Pxx', 1);
  
if nargin < 2
    error(message('signal:bandpower:FreqVectorMissing'));
end

validateattributes(W,{'numeric'},{'vector','finite','real'}, ...
    'bandpower', 'F', 2);
if size(Pxx,1) ~= numel(W)
    error(message('signal:bandpower:FreqVectorMismatch'));
end

if any(diff(W)<0),
    error(message('signal:bandpower:FreqVectorMustBeStrictlyIncreasing'));
end

% Cast to enforce Precision rules
W = double(W);
if nargin < 3,        
    freqrange = [W(1) W(end)];
    freqrangespecified = false;
else
    validateattributes(freqrange,{'numeric'},{'vector','finite','real'}, ...
        'bandpower', 'FREQRANGE', 3);    
    if length(freqrange)~=2
        error(message('signal:bandpower:FreqRangeTwoElement'));
    elseif freqrange(1)<W(1) || freqrange(2)>W(end)
        error(message('signal:bandpower:FreqRangeOutOfBounds'));
    elseif freqrange(1) >= freqrange(2)
        error(message('signal:bandpower:FreqRangeMustBeStrictlyIncreasing'));
    end
    freqrangespecified = true;
end

% force column vector
W = W(:);

% Find indices of freq range requested.
idx1 = find(W<=freqrange(1), 1, 'last' );
idx2 = find(W>=freqrange(2), 1, 'first');

% Determine the width of the rectangle used to approximate the integral.
width = diff(W);
if freqrangespecified
    lastRectWidth = 0;  % Don't include last point of PSD data.
    width = [width; lastRectWidth];
else
    % There are two cases when spectrum is twosided, CenterDC or not.
    % In both cases, the frequency samples does not cover the entire
    % 2*pi (or Fs) region due to the periodicity.  Therefore, the
    % missing freq range has to be compensated in the integral.  The
    % missing freq range can be calculated as the difference between
    % 2*pi (or Fs) and the actual frequency vector span.  For example,
    % considering 1024 points over 2*pi, then frequency vector will be
    % [0 2*pi*(1-1/1024)], i.e., the missing freq range is 2*pi/1024.
    %
    % When CenterDC is true, if the number of points is even, the
    % Nyquist point (Fs/2) is exact, therefore, the missing range is at
    % the left side, i.e., the beginning of the vector.  If the number
    % of points is odd, then the missing freq range is at both ends.
    % However, due to the symmetry of the real signal spectrum, it can
    % still be considered as if it is missing at the beginning of the
    % vector.  Even when the spectrum is asymmetric, since the
    % approximation of the integral is close when NFFT is large,
    % putting it in the beginning of the vector is still ok.
    %
    % When CenterDC is false, the missing range is always at the end of
    % the frequency vector since the frequency always starts at 0.
    
    % assuming a relatively uniform interval
    missingWidth = (W(end) - W(1)) / (numel(W) - 1);
    
    % if CenterDC was not specified, the first frequency point will
    % be 0 (DC).
    centerDC = ~isequal(W(1),0);
    if centerDC
        width = [missingWidth; width];
    else
        width = [width; missingWidth];
    end
end

% Sum the average power over the range of interest.
pwr = width(idx1:idx2)'*Pxx(idx1:idx2,:);

% [EOF]