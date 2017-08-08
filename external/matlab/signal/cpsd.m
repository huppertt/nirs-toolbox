function varargout = cpsd(x,y,varargin)
%CPSD   Cross Power Spectral Density (CPSD) estimate via Welch's method.
%   Pxy = CPSD(X,Y) returns the Cross Power Spectral Density estimate, Pxy,
%   of two discrete-time signals, X and Y, using Welch's averaged,
%   modified periodogram method. By default, X and Y are divided into
%   eight sections with 50% overlap, each section is windowed with a
%   Hamming window and eight modified periodograms are computed and
%   averaged. See "help pwelch" and "help cpsd" for complete details.
%
%   X and Y may be either vectors or two-dimensional matrices. If both are
%   matrices, they must have the same size, and CPSD operates columnwise:
%   Pxy(:,n) = CPSD(X(:,n),Y(:,n)). If one is a matrix and the other is a
%   vector, the vector is converted to a column vector and internally
%   expanded so both inputs have the same number of columns.
%
%   Pxy is the distribution of power per unit frequency. For real signals,
%   CPSD returns the one-sided Cross PSD by default; for complex signals,
%   it returns the two-sided Cross PSD. Note that a one-sided Cross PSD
%   contains the total power of the input signal.
%
%   Pxy = CPSD(X,Y,WINDOW), when WINDOW is a vector, divides each column of
%   X and Y into overlapping sections of length equal to the length of
%   WINDOW, and then windows each section with the vector specified in
%   WINDOW. If WINDOW is an integer, then each column of X and Y are
%   divided into sections of length WINDOW, and each section is windowed
%   with a Hamming of that length. If WINDOW is omitted or specified as
%   empty, a Hamming window is used to obtain eight sections of X and Y.
%
%   Pxy = CPSD(X,Y,WINDOW,NOVERLAP) uses NOVERLAP samples of overlap from
%   section to section. NOVERLAP must be an integer smaller than the
%   length of WINDOW if WINDOW is a vector, or smaller than WINDOW if
%   WINDOW is an integer. If NOVERLAP is omitted or specified as empty, it
%   is set to obtain a 50% overlap.
%
%   [Pxy,W] = CPSD(X,Y,WINDOW,NOVERLAP,NFFT) specifies the number of FFT
%   points used to calculate the Cross PSD estimate. For real signals, Pxy
%   has length (NFFT/2+1) if NFFT is even, and (NFFT+1)/2 if NFFT is odd.
%   For complex signals, Pxy always has length NFFT. If NFFT is specified
%   as empty, the  default NFFT -the maximum of 256 or the next power of
%   two greater than the length of each section of X (and Y)- is used.
%
%   If NFFT is greater than the length of each section, the data is
%   zero-padded. If NFFT is less than the section length, the segment is
%   "wrapped" (using DATAWRAP) to make the length equal to NFFT. This
%   produces the correct FFT when NFFT is smaller than the section length.
%
%   W is the vector of normalized frequencies at which the PSD is
%   estimated. W has units of radians/sample. For real signals, W spans
%   the interval [0,pi] when NFFT is even and [0,pi) when NFFT is odd. For
%   complex signals, W always spans the interval [0,2*pi).
%
%   [Pxy,W] = CPSD(X,Y,WINDOW,NOVERLAP,W) computes the two-sided Cross CPSD
%   at the normalized angular frequencies contained in the vector W. W must
%   have at least two elements.
%
%   [Pxy,F] = CPSD(X,Y,WINDOW,NOVERLAP,NFFT,Fs) returns the Cross PSD as a
%   function of physical frequency. Fs is the sampling frequency specified
%   in hertz. If Fs is empty, it defaults to 1 Hz.
%
%   F is the vector of frequencies (in hertz) at which Pxy is estimated.
%   For real signals, F spans the interval [0,Fs/2] when NFFT is even and
%   [0,Fs/2) when NFFT is odd. For complex signals, F always spans the
%   interval [0,Fs).
%
%   [Pxy,F] = CPSD(X,Y,WINDOW,NOVERLAP,F,Fs) computes the Cross PSD
%   estimate at the physical frequencies contained in the vector F. F must
%   be expressed in hertz and have at least two elements. 
%
%   [...] = CPSD(...,FREQRANGE) returns the Cross PSD computed over the
%   specified range of frequencies based upon the value of FREQRANGE:
%
%      'onesided' - returns the one-sided Cross PSD of real input signals X
%         and Y. If NFFT is even, Pxy has length NFFT/2+1 and is computed
%         over the interval [0,pi]. If NFFT is odd, Pxy has length
%         (NFFT+1)/2 and is computed over the interval [0,pi). When Fs is
%         optionally specified, the intervals become [0,Fs/2) and [0,Fs/2]
%         for even and odd NFFT, respectively.
%
%      'twosided' - returns the two-sided Cross PSD for either real or
%         complex input X and Y. Pxy has length NFFT and is computed over
%         the interval [0,2*pi). When Fs is specified, the interval becomes
%         [0,Fs).
%
%      'centered' - returns the centered two-sided Cross PSD for either
%         real or complex X and Y. Pxy has length NFFT and is computed
%         over the interval (-pi, pi] for even NFFT and (-pi, pi) for odd
%         NFFT. When Fs is specified, the intervals become (-Fs/2, Fs/2]
%         and (-Fs/2, Fs/2) for even and odd NFFT, respectively.
%
%      FREQRANGE may be placed in any position in the input argument list
%      after NOVERLAP. The default value of FREQRANGE is 'onesided' when X
%      and Y are both real and 'twosided' when either X or Y is complex.
%
%   CPSD(...) with no output arguments plots the Cross PSD (in decibels per
%   unit frequency) in the current figure window.
%
%   EXAMPLE:
%      Fs = 1000;   t = 0:1/Fs:.296;
%      x = cos(2*pi*t*200)+randn(size(t));  % A cosine of 200Hz plus noise
%      y = cos(2*pi*t*100)+randn(size(t));  % A cosine of 100Hz plus noise
%      cpsd(x,y,[],[],[],Fs,'twosided');    % Uses default window, overlap & NFFT. 
% 
%   See also PWELCH, PERIODOGRAM, PCOV, PMCOV, PBURG, PYULEAR, PEIG, PMTM,
%   PMUSIC.

%   Copyright 1988-2014 The MathWorks, Inc.

%   References:
%     [1] Petre Stoica and Randolph Moses, Introduction To Spectral
%         Analysis, Prentice-Hall, 1997, pg. 15
%     [2] Monson Hayes, Statistical Digital Signal Processing and 
%         Modeling, John Wiley & Sons, 1996.

narginchk(1,7);
nargoutchk(0,3);

esttype = 'cpsd';
% Possible outputs are:
%       Plot
%       Pxx
%       Pxx, freq
[varargout{1:nargout}] = welch({x,y},esttype,varargin{:});

% [EOF]
