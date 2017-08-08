function varargout = psd(this,x,varargin)
%PSD   Power Spectral Density (PSD) estimate.
%
%   PSD is not recommended.
%   Use the following functions instead:
%      <a href="matlab:help pburg">pburg</a>
%      <a href="matlab:help pcov">pcov</a>
%      <a href="matlab:help pmcov">pmcov</a>
%      <a href="matlab:help pyulear">pyulear</a>
%      <a href="matlab:help periodogram">periodogram</a>
%      <a href="matlab:help pmtm">pmtm</a>
%      <a href="matlab:help pwelch">pwelch</a>
%
%    Hpsd = PSD(H,X) returns a DSP data object (<a href="matlab:help dspdata">dspdata</a>) that contains the
%    PSD estimate of the discrete-time signal vector X estimated using the
%    PSD estimator specified in H.
%
%    Valid PSD estimators:
%           <a href="matlab:help spectrum.periodogram">periodogram</a>    <a href="matlab:help spectrum.mcov">mcov</a>
%           <a href="matlab:help spectrum.welch">welch</a>          <a href="matlab:help spectrum.mtm">mtm</a>
%           <a href="matlab:help spectrum.burg">burg</a>           <a href="matlab:help spectrum.yulear">yulear</a>
%           <a href="matlab:help spectrum.cov">cov</a>
%
%    The PSD data contained in the object Hpsd is the distribution of power
%    per unit frequency. For real signals, PSD returns the one-sided PSD by
%    default; for complex signals, it returns the two-sided PSD.  Note that
%    a one-sided PSD contains the total power of the input signal.
%
%    The power spectral density is intended for continuous spectra. Note
%    that unlike the mean-squared spectrum (MSS), in this case the peaks in
%    the spectra do not reflect the power at a given frequency. Instead,
%    the integral of the PSD over a given frequency band computes the
%    average power in the signal over such frequency band. See the help on
%    AVGPOWER for more information.
% 
%    The Hpsd object also contains a vector of normalized frequencies W at
%    which the PSD is estimated.  W has units of rad/sample.  For real
%    signals, W spans the interval [0,Pi] when NFFT is even and [0,Pi) when
%    NFFT is odd.  For complex signals, W always spans the interval
%    [0,2*Pi).
% 
%    Hpsd = PSD(...,'Fs',Fs) returns a PSD object with the spectrum
%    computed as a function of physical frequency (Hz).  Fs is the sampling
%    frequency specified in Hz.
% 
%    Hpsd = PSD(...,'SpectrumType','twosided') returns a PSD object with a
%    two-sided PSD of a real signal X. In this case, the spectrum will be
%    computed over the interval [0,2*Pi) if Fs is not specified and over
%    the interval [0,Fs) if Fs is specified.  The SpectrumType can also be
%    'onesided' for a real signal X, which is the default behavior.
% 
%    Hpsd = PSD(...,'NFFT',nfft) specifies nfft as the number of FFT points
%    to use to calculate the power spectral density.
%
%    Hpsd = PSD(...,'CenterDC',true) specifies that the spectrum should be
%    shifted so that the zero-frequency component is in the center of the
%    spectrum.  CenterDC is set to false by default.
%
%    Hpsd = PSD(...,'FreqPoints','User Defined','FrequencyVector',f) 
%    returns a PSD object evaluated at the frequencies defined by the
%    vector f of frequencies
%     
%    Hpsd = PSD(...,'ConfLevel',p) specifies the confidence level p. p
%    takes values in the interval (0,1). The confidence interval is added
%    to the Hpsd object as a two column matrix. The first column contains
%    the lower bound while the second column contains the upper bound. 
%
%    PSD(...) with no output arguments by default plots the PSD estimate in
%    dB per unit frequency in the current figure window.
%  
%    An alternative to specifying the individual input arguments to PSD is
%    to create an options object using <a href="matlab:help spectrum/psdopts">psdopts</a>.
%
%   EXAMPLE : Spectral analysis of a signal that contains a 200Hz cosine
%             % plus noise.
%             Fs = 1000;   t = 0:1/Fs:.296;
%             x = cos(2*pi*t*200)+randn(size(t));  
%             h = spectrum.welch;                  % Create a Welch spectral estimator. 
%             psd(h,x,'Fs',Fs);                    % Calculate and plot the one-sided PSD.
%             hpsd = psd(h,x,'ConfLevel',.98);     % PSD with confidence level
%             figure,plot(hpsd)

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help file for PSD method.

% [EOF]
