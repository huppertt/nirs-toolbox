function varargout = msspectrum(this,x,varargin)
%MSSPECTRUM   Mean-square spectrum (MSS).
%
%   MSSPECTRUM is not recommended.  Use <a href="matlab:help periodogram">periodogram</a> or <a href="matlab:help pwelch">pwelch</a> instead.
%
%    Hmss = MSSPECTRUM(H,X) returns a DSP data object (<a href="matlab:help dspdata">dspdata</a>) in Hmss
%    that contains the mean-square (power) estimate of the discrete-time
%    signal vector X estimated using the MSS estimator specified in H.
%
%    Valid MSS estimators:
%           <a href="matlab:help spectrum.periodogram">periodogram</a>
%           <a href="matlab:help spectrum.welch">welch</a>
%
%    The mean-square spectrum contained in the object Hmss is the
%    distribution of power over frequency. For real signals, MSSPECTRUM
%    returns the one-sided mean-square spectrum by default; for complex
%    signals, it returns the two-sided mean-square spectrum.  Note that a
%    one-sided mean-square spectrum contains the total power of the input
%    signal.
% 
%    The mean-squared spectrum is intended for discrete spectra. Unlike the
%    power spectral density (PSD), the peaks in the mean-square spectrum
%    reflect the power in the signal at a given frequency.
% 
%    The Hmss object also contains the vector of normalized frequencies W
%    at which the mean-square spectrum is estimated.  W has units of
%    rad/sample.  For real signals, W spans the interval [0,Pi] when NFFT
%    is even and [0,Pi) when NFFT is odd.  For complex signals, W always
%    spans the interval [0,2*Pi).
% 
%    Hmss = MSSPECTRUM(H,X,'Fs',Fs) returns a mean-square spectrum object
%    containing the spectrum computed as a function of physical frequency
%    (Hz). Fs is the sampling frequency specified in Hz.
% 
%    Hmss = MSSPECTRUM(...,'SpectrumType','twosided') returns a mean-square
%    spectrum object with a two-sided mean-square spectrum of a real signal
%    X. In this case, the spectrum will be computed over the interval
%    [0,2*Pi) if Fs is not specified and over the interval [0,Fs) if Fs is
%    specified. The SpectrumType can also be 'onesided' for a real signal
%    X, which is the default behavior.
% 
%    Hmss = MSSPECTRUM(...,'NFFT',nfft) specifies nfft as the number of FFT
%    points to use to calculate the mean-square spectrum.
%
%    Hmss = MSSPECTRUM(...,'CenterDC',true) specifies that the spectrum
%    should be shifted, so that the zero-frequency component is in the
%    center of the spectrum.  The CenterDC property is by default set to
%    false.
%
%    Hmss = MSSPECTRUM(...,'FreqPoints','User Defined','FrequencyVector',f) 
%    returns a mean-square spectrum object evaluated at the frequencies
%    defined by the vector f of frequencies.
%
%    Hmss = MSSPECTRUM(...,'ConfLevel',p) specifies the confidence level p.
%    p takes values in the interval (0,1). The confidence interval is added
%    to the Hmss object as a two column matrix. The first column contains
%    the lower bound while the second column contains the upper bound. 
%
%    MSSPECTRUM(...) with no output arguments by default plots the
%    mean-square spectrum estimate in dB per unit frequency in the current
%    figure window.
%  
%    An alternative to specifying the individual input arguments to
%    MSSPECTRUM is to create an options object using <a href="matlab:help spectrum/msspectrumopts">msspectrumopts</a>.
%
%   EXAMPLE 1: 
%      % In this example we will measure the power of a deterministic
%      % power signal which has a frequency component at 200Hz. We'll use a
%      % signal with a peak amplitude of 3 volts therefore, the theoretical
%      % power at 200Hz should be 3^2/2 volts^2 (watts) or 6.5321dB. 
%      Fs = 1000;   t = 0:1/Fs:.2;
%      x = 3*cos(2*pi*t*200);  
%      h = spectrum.welch;           % Create a Welch spectral estimator. 
%
%      % Plot the one-sided mean-square spectrum.
%      msspectrum(h,x,'Fs',Fs,'NFFT',2^14);
%
%   EXAMPLE 2: Evaluate the mean-square spectrum at frequencies around 200 Hz.
%      msspectrum(h,x,'Fs',Fs,'FreqPoints','User Defined','SpectrumType',...
%      'TwoSided','FrequencyVector',linspace(150,250,100),'ConfLevel',0.95);

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help file for MSSPECTRUM method.

% [EOF]
