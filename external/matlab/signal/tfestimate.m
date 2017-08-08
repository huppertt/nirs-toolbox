function varargout = tfestimate(x,y,varargin)
%TFESTIMATE   Transfer Function Estimate.
%   Txy = TFESTIMATE(X,Y) estimates the transfer function of the system
%   with input X and output Y using Welch's averaged, modified periodogram
%   method. Txy is the quotient of the Cross Power Spectral Density (CPSD)
%   of X and Y, Pxy, and the Power Spectral Density (PSD) of X, Pxx.
%
%   By default, X and Y are divided into eight sections with 50% overlap,
%   each section is windowed with a Hamming window, and eight modified
%   periodograms are computed and averaged.  See "help pwelch" and "help
%   cpsd" for complete details.
%
%   X and Y may be either vectors or two-dimensional matrices.  If both are
%   matrices, they must have the same size and TFESTIMATE operates
%   columnwise: Txy(:,n) = TFESTIMATE(X(:,n),Y(:,n)).  If one is a matrix
%   and the other is a vector, the vector is converted to a column vector
%   and internally expanded so both inputs have the same number of columns.
%
%   Txy = TFESTIMATE(X,Y,WINDOW), when WINDOW is a vector, divides each
%   column of X and Y into overlapping sections of length equal to the
%   length of WINDOW, and then windows each section with the vector
%   specified in WINDOW.  If WINDOW is an integer, then each column of X
%   and Y are divided into sections of length WINDOW, and each section is
%   windowed with a Hamming of that length.  If WINDOW is omitted or
%   specified as empty, a Hamming window is used to obtain eight sections
%   of X and Y.
%    
%   Txy = TFESTIMATE(X,Y,WINDOW,NOVERLAP) uses NOVERLAP samples of overlap
%   from section to section.  NOVERLAP must be an integer smaller than the
%   WINDOW if WINDOW is an integer, or smaller than the length of WINDOW if
%   WINDOW is a vector.  If NOVERLAP is omitted or specified as empty, it
%   is set to obtain a 50% overlap.
%
%   [Txy,W] = TFESTIMATE(X,Y,WINDOW,NOVERLAP,NFFT) specifies the number of
%   FFT points used to calculate the PSD and CPSD estimates.  For real X
%   and Y, Txy has length (NFFT/2+1) if NFFT is even, and (NFFT+1)/2 if
%   NFFT is odd.  For complex X or Y, Txy always has length NFFT.  If NFFT
%   is specified as empty, it is set to either 256 or the next power of two
%   greater than the length of each section of X (or Y), whichever is
%   greater.
%
%   If NFFT is greater than the length of each section, the data is
%   zero-padded. If NFFT is less than the section length, the segment is
%   "wrapped" (using DATAWRAP) to make the length equal to NFFT. This
%   produces the correct FFT when NFFT is smaller than the section length.
%
%   W is the vector of normalized frequencies at which Txy is estimated.
%   W has units of radians/sample.  For real signals, W spans the interval
%   [0,pi] when NFFT is even and [0,pi) when NFFT is odd.  For complex
%   signals, W always spans the interval [0,2*pi).
%
%   [Txy,W] = TFESTIMATE(X,Y,WINDOW,NOVERLAP,W) computes the two-sided 
%   transfer function estimate at the normalized angular frequencies 
%   contained in the vector W.  W must have at least two elements.
%
%   [Txy,F] = TFESTIMATE(X,Y,WINDOW,NOVERLAP,NFFT,Fs) returns the transfer
%   function estimate as a function of physical frequency.  Fs is the
%   sampling frequency specified in hertz.  If Fs is empty, it defaults to
%   1 Hz.
%
%   F is the vector of frequencies (in hertz) at which Txy is estimated.
%   For real signals, F spans the interval [0,Fs/2] when NFFT is even and
%   [0,Fs/2) when NFFT is odd.  For complex signals, F always spans the
%   interval [0,Fs).
%
%   [Txy,F] = TFESTIMATE(X,Y,WINDOW,NOVERLAP,F,Fs) computes the transfer
%   function estimate at the physical frequencies contained in the vector
%   F.  F must be expressed in hertz and have at least two elements.
%
%   [...] = TFESTIMATE(...,FREQRANGE) returns the transfer function
%   estimate computed over the specified range of frequencies based upon
%   the value of FREQRANGE:
%
%      'onesided' - returns the one-sided transfer function estimate of
%         real input signals X and Y. If NFFT is even, Txy has length
%         NFFT/2+1 and is computed over the interval [0,pi]. If NFFT is
%         odd, Txy has length (NFFT+1)/2 and is computed over the interval
%         [0,pi). When Fs is optionally specified, the intervals become
%         [0,Fs/2) and [0,Fs/2] for even and odd NFFT, respectively.
%
%      'twosided' - returns the two-sided transfer function estimate for
%         either real or complex input X and Y. Txy has length NFFT and is
%         computed over the interval [0,2*pi). When Fs is specified, the
%         interval becomes [0,Fs).
%
%      'centered' - returns the centered two-sided two-sided transfer
%         function estimate for either real or complex X and Y.  Txy has
%         length NFFT and is computed over the interval (-pi, pi] for even
%         NFFT and (-pi, pi) for odd NFFT. When Fs is specified, the
%         intervals become (-Fs/2, Fs/2] and (-Fs/2, Fs/2) for even and odd
%         NFFT, respectively.
%
%      FREQRANGE may be placed in any position in the input argument list
%      after NOVERLAP.  The default value of FREQRANGE is 'onesided' when X
%      and Y are both real and 'twosided' when either X or Y is complex.
%
%   TFESTIMATE(...) with no output arguments plots the transfer function
%   estimate (in decibels per unit frequency) in the current figure window.
%
%   EXAMPLE:
%      h = fir1(30,0.2,rectwin(31));
%      x = randn(16384,1);
%      y = filter(h,1,x);
%      tfestimate(x,y,[],512,1024); % Plot estimate using a default window.
% 
%      % Estimate the transfer function and calculate the
%      % filter coefficients using INVFREQZ.
%      [txy,w] = tfestimate(x,y,[],512,1024); 
%      [b,a] = invfreqz(txy,w,30,0);
%      htool = fvtool(h,1,b,a); legend(htool,'Original','Estimated');
%
%   See also CPSD, PWELCH, MSCOHERE, PERIODOGRAM. 

%   Copyright 1988-2014 The MathWorks, Inc.

narginchk(2,7)

esttype = 'tfe';
% Possible outputs are:
%       Plot
%       Txy
%       Txy, freq
[varargout{1:nargout}] = welch({x,y},esttype,varargin{:});

if nargout == 0,
    title('Transfer Function Estimate via Welch');
end

% [EOF]
