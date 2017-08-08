function [Px,w,Pxxc] = periodogram(x,varargin)
%PERIODOGRAM  Power Spectral Density (PSD) estimate via periodogram method.
%   Pxx = PERIODOGRAM(X) returns the PSD estimate, Pxx, of a signal, X.
%   When X is a vector, it is converted to a column vector and treated as a
%   single channel.  When X is a matrix, the PSD is computed independently
%   for each column and stored in the corresponding column of Pxx.
%  
%   By default, the signal X is windowed with a rectangular window of the
%   same length as X. The PSD estimate is computed using an FFT of length
%   given by the larger of 256 and the next power of 2 greater than the
%   length of X.
%
%   Note that the default window (rectangular) has a 13.3 dB sidelobe
%   attenuation. This may mask spectral content below this value (relative
%   to the peak spectral content). Choosing different windows will enable
%   you to make tradeoffs between resolution (e.g., using a rectangular
%   window) and sidelobe attenuation (e.g., using a Hann window). See
%   WinTool for more details.
%
%   Pxx is the distribution of power per unit frequency. For real signals,
%   PERIODOGRAM returns the one-sided PSD by default; for complex signals,
%   it returns the two-sided PSD.  Note that a one-sided PSD contains the
%   total power of the input signal.
%
%   Pxx = PERIODOGRAM(X,WINDOW) specifies a window to be applied to X. If X
%   is a vector, WINDOW must be a vector of the same length as X.  If X is
%   a matrix, WINDOW must be a vector whose length is the same as the
%   number of rows of X.  If WINDOW is a window other than a rectangular,
%   the resulting estimate is a modified periodogram. If WINDOW is
%   specified as empty, the default window is used.
%
%   Pxx = PERIODOGRAM(X,WINDOW,...,SPECTRUMTYPE) uses the window scaling
%   algorithm specified by SPECTRUMTYPE when computing the power spectrum:
%     'psd'   - returns the power spectral density
%     'power' - scales each estimate of the PSD by the equivalent noise
%               bandwidth (in Hz) of the window.  Use this option to
%               obtain an estimate of the power at each frequency.
%   The default value for SPECTRUMTYPE is 'psd'
%
%   [Pxx,W] = PERIODOGRAM(X,WINDOW,NFFT) specifies the number of FFT points
%   used to calculate the PSD estimate.  For real X, Pxx has length
%   (NFFT/2+1) if NFFT is even, and (NFFT+1)/2 if NFFT is odd.  For complex
%   X, Pxx always has length NFFT.  If NFFT is specified as empty, the
%   default NFFT is used.
%
%   Note that if NFFT is greater than the length of WINDOW, the data is
%   zero-padded. If NFFT is less than the length of WINDOW, the segment is
%   "wrapped" (using DATAWRAP) to make the length equal to NFFT to produce
%   the correct FFT.
%
%   W is the vector of normalized frequencies at which the PSD is
%   estimated.  W has units of radians/sample.  For real signals, W spans
%   the interval [0,pi] when NFFT is even and [0,pi) when NFFT is odd.  For
%   complex signals, W always spans the interval [0,2*pi).
%
%   [Pxx,W] = PERIODOGRAM(X,WINDOW,W) computes the two-sided PSD at the
%   normalized angular frequencies contained in the vector W. W must have
%   at least two elements.
%
%   [Pxx,F] = PERIODOGRAM(X,WINDOW,NFFT,Fs) returns a PSD computed as
%   a function of physical frequency.  Fs is the sampling frequency
%   specified in hertz.  If Fs is empty, it defaults to 1 Hz.
%
%   F is the vector of frequencies (in hertz) at which the PSD is
%   estimated.  For real signals, F spans the interval [0,Fs/2] when NFFT
%   is even and [0,Fs/2) when NFFT is odd.  For complex signals, F always
%   spans the interval [0,Fs).
%
%   [Pxx,F] = PERIODOGRAM(X,WINDOW,F,Fs) computes the two-sided PSD at the 
%   frequencies contained in vector F. F must contain at least two elements
%   and be expressed in units of hertz.
%
%   [...] = PERIODOGRAM(X,WINDOW,NFFT,...,FREQRANGE) returns the PSD
%   over the specified range of frequencies based upon the value of
%   FREQRANGE:
%
%      'onesided' - returns the one-sided PSD of a real input signal X.
%         If NFFT is even, Pxx has length NFFT/2+1 and is computed over the
%         interval [0,pi].  If NFFT is odd, Pxx has length (NFFT+1)/2 and
%         is computed over the interval [0,pi). When Fs is specified, the
%         intervals become [0,Fs/2) and [0,Fs/2] for even and odd NFFT,
%         respectively.
%
%      'twosided' - returns the two-sided PSD for either real or complex
%         input X.  Pxx has length NFFT and is computed over the interval
%         [0,2*pi). When Fs is specified, the interval becomes [0,Fs).
%
%      'centered' - returns the centered two-sided PSD for either real or
%         complex X.  Pxx has length NFFT and is computed over the interval
%         (-pi, pi] for even length NFFT and (-pi, pi) for odd length NFFT.
%         When Fs is specified, the intervals become (-Fs/2, Fs/2] and
%         (-Fs/2, Fs/2) for even and odd NFFT, respectively.
%
%      FREQRANGE may be placed in any position in the input argument list
%      after WINDOW.  The default value of FREQRANGE is 'onesided' when X
%      is real and 'twosided' when X is complex.
%
%   [Pxx,F,Pxxc] = PERIODOGRAM(...,'ConfidenceLevel',P) , where P is a
%   scalar between 0 and 1, returns the P*100% confidence interval for Pxx.
%   The default value for P is .95.  Confidence intervals are computed
%   using a chi-squared approach.  Pxxc has twice as many columns as Pxx.
%   Odd-numbered columns contain the lower bounds of the confidence
%   intervals; even-numbered columns contain the upper bounds.  Thus,
%   Pxxc(M,2*N-1) is the lower bound and Pxxc(M,2*N) is the upper bound
%   corresponding to the estimate Pxx(M,N).
%
%   PERIODOGRAM(...) with no output arguments by default plots the PSD
%   estimate (in decibels per unit frequency) in the current figure window.
%
%   EXAMPLE:
%      % Compute the two-sided periodogram of a 200 Hz sinusoid embedded
%      % in noise using the hamming window.
%      Fs = 1000;   t = 0:1/Fs:.3;
%      x = cos(2*pi*t*200)+randn(size(t));
%      periodogram(x,[],'twosided',512,Fs);
%
%   See also PWELCH, PBURG, PCOV, PYULEAR, PMTM, PMUSIC, PMCOV, PEIG.

%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(1,9);

% look for psd, power, and ms window compensation flags
[esttype, varargin] = psdesttype({'psd','power','ms'},'psd',varargin);

if isvector(x)
    N = length(x); % Record the length of the data
else
    N = size(x,1);
end

% extract window argument
if ~isempty(varargin) && ~ischar(varargin{1})
    win = varargin{1};
    varargin = varargin(2:end);
else
    win = [];
end

% Generate a default window if needed
winName = 'User Defined';
winParam = '';
if isempty(win),
    win = rectwin(N);
    winName = 'Rectangular';
    winParam = N;
end

% Cast to enforce precision rules
if any([signal.internal.sigcheckfloattype(x,'single','periodogram','X')...
    signal.internal.sigcheckfloattype(win,'single','periodogram','WINDOW')]) 
  x = single(x);
  win = single(win);
end

[options,msg,msgobj] = periodogram_options(isreal(x),N,varargin{:});
if ~isempty(msg)
  error(msgobj)
end

Fs    = options.Fs;
nfft  = options.nfft;

% Compute the PS using periodogram over the whole nyquist range.
[Sxx,w] = computeperiodogram(x,win,nfft,esttype,Fs);

nrow = 1;
% If frequency vector was specified, return and plot two-sided PSD
% The computepsd function expects NFFT to be a scalar
if (length(nfft) > 1),
    [ncol,nrow] = size(nfft);
    nfft = max(ncol,nrow);
    if (length(options.nfft)>1 && strcmpi(options.range,'onesided'))
        warning(message('signal:periodogram:InconsistentRangeOption'));
        options.range = 'twosided';
    end
end

% Compute the 1-sided or 2-sided PSD [Power/freq] or mean-square [Power].
% Also, compute the corresponding freq vector & freq units.
[Pxx,w,units] = computepsd(Sxx,w,options.range,nfft,Fs,esttype);

% compute confidence intervals if needed.
if ~strcmp(options.conflevel,'omitted')
    Pxxc = confInterval(options.conflevel, Pxx, x, w, options.Fs);
elseif nargout>2
    Pxxc = confInterval(0.95, Pxx, x, w, options.Fs);
else
    Pxxc = [];
end

if nargout==0, % Plot when no output arguments are specified
    w = {w};
    if strcmpi(units,'Hz'), w = [w, {'Fs',options.Fs}]; end
    
    if strcmp(esttype,'psd')
        hdspdata = dspdata.psd(Pxx,w{:},'SpectrumType',options.range);
    else
        hdspdata = dspdata.msspectrum(Pxx,w{:},'SpectrumType',options.range);
    end
    % plot the confidence levels if conflevel is specified.
    if ~isempty(Pxxc)
        hdspdata.ConfLevel = options.conflevel;
        hdspdata.ConfInterval = Pxxc;
    end
    % Create a spectrum object to store in the PSD object's metadata.
    hspec = spectrum.periodogram({winName,winParam});
    hdspdata.Metadata.setsourcespectrum(hspec);
    
    if options.centerdc
        centerdc(hdspdata);
    end
    plot(hdspdata);
    
    if strcmp(esttype,'power')
        title(getString(message('signal:periodogram:PeriodogramPowerSpectrumEstimate')));
    end
else
    if options.centerdc
        [Pxx, w, Pxxc] = psdcenterdc(Pxx, w, Pxxc, options);
    end
    Px = Pxx;
    
    % If the input is a vector and a row frequency vector was specified,
    % return output as a row vector for backwards compatibility.
    if nrow > 1 && isvector(x)
        Px = Px.'; w = w.';% Sxx = Sxx.';
    end
    
    % Cast to enforce precision rules
    % Only case if output is requested, otherwise plot using double
    % precision frequency vector.
    if isa(Px,'single')
      w = single(w);
    end
end

%------------------------------------------------------------------------------
function [options,msg,msgobj] = periodogram_options(isreal_x,N,varargin)
%PERIODOGRAM_OPTIONS   Parse the optional inputs to the PERIODOGRAM function.
%   PERIODOGRAM_OPTIONS returns a structure, OPTIONS, with following fields:
%
%   options.nfft         - number of freq. points at which the psd is estimated
%   options.Fs           - sampling freq. if any
%   options.range        - 'onesided' or 'twosided' psd
%   options.centerdc     - true if 'centered' specified

% Generate defaults
options.nfft = max(256, 2^nextpow2(N));
options.Fs = []; % Work in rad/sample

% Determine if frequency vector specified
freqVecSpec = false;
if (~isempty(varargin) && length(varargin{1}) > 1)
    freqVecSpec = true;
end

if isreal_x && ~freqVecSpec,
    options.range = 'onesided';
else
    options.range = 'twosided';
end

if any(strcmp(varargin, 'whole'))
    warning(message('signal:periodogram:invalidRange', 'whole', 'twosided'));
elseif any(strcmp(varargin, 'half'))
    warning(message('signal:periodogram:invalidRange', 'half', 'onesided'));
end

[options,msg,msgobj] = psdoptions(isreal_x,options,varargin{:});

% Cast to enforce precision rules
options.Fs = double(options.Fs);
options.nfft = double(options.nfft);

%--------------------------------------------------------------------------
function Pxxc = confInterval(CL, Pxx, x, w, fs)
%   Reference: D.G. Manolakis, V.K. Ingle and S.M. Kagon,
%   Statistical and Adaptive Signal Processing,
%   McGraw-Hill, 2000, Chapter 5

% Compute confidence intervals using double precision arithmetic
Pxx = double(Pxx);
x = double(x);

k = 1;
c = chi2conf(CL,k);
PxxcLower = Pxx*c(1);
PxxcUpper = Pxx*c(2);
Pxxc = reshape([PxxcLower; PxxcUpper],size(Pxx,1),2*size(Pxx,2));

% DC and Nyquist bins have only one degree of freedom for real signals
if isreal(x)
    realConf = chi2conf(CL,k/2);
    Pxxc(w == 0,1:2:end) = Pxx(w == 0,:) * realConf(1);
    Pxxc(w == 0,2:2:end) = Pxx(w == 0,:) * realConf(2);
    if isempty(fs)
        Pxxc(w == pi,1:2:end) = Pxx(w == pi,:) * realConf(1);
        Pxxc(w == pi,2:2:end) = Pxx(w == pi,:) * realConf(2);
    else
        Pxxc(w == fs/2,1:2:end) = Pxx(w == fs/2,:) * realConf(1);
        Pxxc(w == fs/2,2:2:end) = Pxx(w == fs/2,:) * realConf(2);
    end
end

% [EOF] periodogram.m
