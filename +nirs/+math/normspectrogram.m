function [Sall,f,t,Pall,fcorr,tcorr] = normspectrogram(x,varargin)
%SPECTROGRAM Spectrogram using a Short-Time Fourier Transform (STFT).
%   S = SPECTROGRAM(X) returns the spectrogram of the signal specified by
%   vector X in the matrix S. By default, X is divided into eight segments
%   with 50% overlap, each segment is windowed with a Hamming window. The
%   number of frequency points used to calculate the discrete Fourier
%   transforms is equal to the maximum of 256 or the next power of two
%   greater than the length of each segment of X.
%
%   If X cannot be divided exactly into eight segments, X will be truncated
%   accordingly.
%
%   S = SPECTROGRAM(X,WINDOW) when WINDOW is a vector, divides X into
%   segments of length equal to the length of WINDOW, and then windows each
%   segment with the vector specified in WINDOW.  If WINDOW is an integer,
%   X is divided into segments of length equal to that integer value, and a
%   Hamming window of equal length is used.  If WINDOW is not specified, the
%   default is used.
%
%   S = SPECTROGRAM(X,WINDOW,NOVERLAP) NOVERLAP is the number of samples
%   each segment of X overlaps. NOVERLAP must be an integer smaller than
%   WINDOW if WINDOW is an integer.  NOVERLAP must be an integer smaller
%   than the length of WINDOW if WINDOW is a vector.  If NOVERLAP is not
%   specified, the default value is used to obtain a 50% overlap.
%
%   S = SPECTROGRAM(X,WINDOW,NOVERLAP,NFFT) specifies the number of
%   frequency points used to calculate the discrete Fourier transforms.
%   If NFFT is not specified, the default NFFT is used.
%
%   S = SPECTROGRAM(X,WINDOW,NOVERLAP,NFFT,Fs) Fs is the sampling frequency
%   specified in Hz. If Fs is specified as empty, it defaults to 1 Hz. If
%   it is not specified, normalized frequency is used.
%
%   Each column of S contains an estimate of the short-term, time-localized
%   frequency content of the signal X.  Time increases across the columns
%   of S, from left to right.  Frequency increases down the rows, starting
%   at 0.  If X is a length NX complex signal, S is a complex matrix with
%   NFFT rows and k = fix((NX-NOVERLAP)/(length(WINDOW)-NOVERLAP)) columns.
%   For real X, S has (NFFT/2+1) rows if NFFT is even, and (NFFT+1)/2 rows
%   if NFFT is odd.
%
%   [S,F,T] = SPECTROGRAM(...) returns a vector of frequencies F and a
%   vector of times T at which the spectrogram is computed. F has length
%   equal to the number of rows of S. T has length k (defined above) and
%   its value corresponds to the center of each segment.
%
%   [S,F,T] = SPECTROGRAM(X,WINDOW,NOVERLAP,F,Fs) computes the two-sided
%   spectrogram at the frequencies specified in vector F.  F must be
%   expressed in hertz and have at least two elements.
%
%   [S,F,T,P] = SPECTROGRAM(...) P is a matrix representing the Power
%   Spectral Density (PSD) of each segment. For real signals, SPECTROGRAM
%   returns the one-sided modified periodogram estimate of the PSD of each
%   segment; for complex signals and in the case when a vector of
%   frequencies is specified, it returns the two-sided PSD.
%
%   [S,F,T,P] = SPECTROGRAM(...,'MinThreshold',THRESH) sets the elements of
%   P to zero when the corresponding elements of 10*log10(P) are less than
%   THRESH. Specify THRESH in units of decibels. The default value of
%   THRESH is -Inf.
%
%   [S,F,T,P] = SPECTROGRAM(...,'reassigned') reassigns each PSD estimate
%   to the location of its center of gravity.  The reassignment is done
%   in-place and returned in P.
%
%   [S,F,T,P,Fc,Tc] = SPECTROGRAM(...) returns the locations in frequency
%   and time of the center of gravity of each estimate in the spectrogram.
%   The frequencies and times are returned in matrices, Fc, and Tc,
%   respectively.  Fc and Tc have the same dimensions as the spectrogram, S.
%
%   [...]  = SPECTROGRAM(...,SPECTRUMTYPE) uses the window scaling
%   algorithm specified by SPECTRUMTYPE when computing the power spectral
%   density matrix P.
%   SPECTRUMTYPE can be set to 'psd' or 'power':
%     'psd'   - returns the power spectral density
%     'power' - scales each estimate of the PSD by the equivalent noise
%               bandwidth of the window (in hertz).  Use this option to
%               obtain an estimate of the power at each frequency.
%   The default value for SPECTRUMTYPE is 'psd'.
%
%   [...] = SPECTROGRAM(...,FREQRANGE)  returns the PSD over the specified
%   range of frequencies based upon the value of FREQRANGE:
%
%      'onesided' - returns the one-sided matrix P of a real input signal X.
%         If NFFT is even, P has length NFFT/2+1 and is computed over the
%         interval [0,pi].  If NFFT is odd, P has length (NFFT+1)/2 and
%         is computed over the interval [0,pi). When Fs is specified, the
%         intervals become [0,Fs/2) and [0,Fs/2] for even and odd NFFT,
%         respectively.
%
%      'twosided' - returns the two-sided matrix P for either real or complex
%         input X.  P has length NFFT and is computed over the interval
%         [0,2*pi). When Fs is specified, the interval becomes [0,Fs).
%
%      'centered' - returns the centered two-sided matrix P for either real
%         or complex X.  P has length NFFT and is computed over the interval
%         (-pi, pi] for even length NFFT and (-pi, pi) for odd length NFFT.
%         When Fs is specified, the intervals become (-Fs/2, Fs/2] and
%         (-Fs/2, Fs/2) for even and odd NFFT, respectively.
%
%      FREQRANGE may be placed in any position in the input argument list
%      after NOVERLAP.  The default value of FREQRANGE is 'onesided' when X
%      is real and 'twosided' when X is complex.
%
%   SPECTROGRAM(...) with no output arguments plots the PSD estimate for
%   each segment on a surface in the current figure. It uses
%   SURF(f,t,10*log10(abs(P)) where P is the fourth output argument.
%
%   SPECTROGRAM(...,FREQLOCATION) controls where MATLAB displays the
%   frequency axis on the plot. This string can be either 'xaxis' or
%   'yaxis'.  Setting this FREQLOCATION to 'yaxis' displays frequency on
%   the y-axis and time on the x-axis.  The default is 'xaxis' which
%   displays the frequency on the x-axis. If FREQLOCATION is specified when
%   output arguments are requested, it is ignored.
%
%   EXAMPLE 1: Spectrogram of quadratic chirp
%     t=0:0.001:2;                    % 2 secs @ 1kHz sample rate
%     y=chirp(t,100,1,200,'q');       % Start @ 100Hz, cross 200Hz at t=1sec
%     spectrogram(y,kaiser(128,18),120,128,1E3,'yaxis');
%     title('Quadratic Chirp: start at 100Hz and cross 200Hz at t=1sec');
%
%   EXAMPLE 2: Reassigned spectrogram of quadratic chirp
%     t=0:0.001:2;                    % 2 secs @ 1kHz sample rate
%     y=chirp(t,100,1,200,'q');       % Start @ 100Hz, cross 200Hz at t=1sec
%     spectrogram(y,kaiser(128,18),120,128,1E3,'reassigned','yaxis');
%     title('Quadratic Chirp: start at 100Hz and cross 200Hz at t=1sec');
%
%   EXAMPLE 3:  Plot instantaneous frequency of quadratic chirp
%     t=0:0.001:2;                    % 2 secs @ 1kHz sample rate
%     y=chirp(t,100,1,200,'q');       % Start @ 100Hz, cross 200Hz at t=1sec
%     % remove estimates less than -30 dB
%     [~,~,~,P,Fc,Tc] = spectrogram(y,kaiser(128,18),120,128,1E3,'minthreshold',-30);
%     plot(Tc(P>0),Fc(P>0),'. ')
%     title('Quadratic Chirp: start at 100Hz and cross 200Hz at t=1sec');
%     xlabel('Time (s)')
%     ylabel('Frequency (Hz)')
%
%   EXAMPLE 4: Waterfall display of the PSD of each segment of a VCO
%     Fs = 10e3;
%     t = 0:1/Fs:2;
%     x1 = vco(sawtooth(2*pi*t,0.5),[0.1 0.4]*Fs,Fs);
%     spectrogram(x1,kaiser(256,5),220,512,Fs);
%     view(-45,65)
%     colormap bone

%   See also PERIODOGRAM, PWELCH, SPECTRUM, GOERTZEL.

% [1] Oppenheim, A.V., and R.W. Schafer, Discrete-Time Signal Processing,
%     Prentice-Hall, Englewood Cliffs, NJ, 1989, pp. 713-718.
% [2] Mitra, S. K., Digital Signal Processing. A Computer-Based Approach.
%     2nd Ed. McGraw-Hill, N.Y., 2001.
% [3] Chassande-Mottin, E., Auger F., and Flandrin, P., Reassignment.
%     Chapter 9.3.1 in Time-Frequency Analysis: Concepts and Methods. F.
%     Hlawatsch and F. Auger Eds. John Wiley & Sons, 2008, pg. 258-259

%   Copyright 1988-2015 The MathWorks, Inc.


maxcomp=30;
n=round(size(x,1)/30);
[X]=convmtx(x,n);
[U,S,V]=nirs.math.mysvd(X(n+1:end-n,:));

nn=min(find(cumsum(diag(S))/sum(diag(S))>.99));
nn=min(nn,maxcomp);

[icasig, A, W] = fastica(X(n+1:end-n,:)', 'numOfIC',nn);

i=1;
[s,f,t,Pxx,fcorr,tcorr] = spectrogram(icasig(i,:),varargin{:});
s=s/sqrt(sum(s(:).*conj(s(:))))/length(t);
Pxx=Pxx/sqrt(sum(Pxx(:).*conj(Pxx(:))))/length(t);
Sall=s/size(icasig,1);
Pall=Pxx/size(icasig,1);

for i=2:size(icasig,1); 
    [s,f,t,Pxx,fcorr,tcorr] = spectrogram(icasig(i,:),varargin{:});
    s=s/sqrt(sum(s(:).*conj(s(:))))/length(t);
    Pxx=Pxx/sqrt(sum(Pxx(:).*conj(Pxx(:))))/length(t);
    Sall=Sall+s/size(icasig,1);
    Pall=Pall+Pxx/size(icasig,1);
end


[s,f,t,Pxx,fcorr,tcorr] = spectrogram(x,varargin{:});

Sall=Sall/sqrt(sum(Sall(:).*conj(Sall(:))))*sqrt(sum(s(:).*conj(s(:))));
Pall=Pall/sqrt(sum(Pall(:).*conj(Pall(:))))*sqrt(sum(Pxx(:).*conj(Pxx(:))));




