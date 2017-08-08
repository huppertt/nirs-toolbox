function hopts = psdopts(this,x)
%PSDOPTS   Power spectral density (PSD) options object.
%
%   PSDOPTS is not recommended.
%   Use the following functions instead:
%      <a href="matlab:help pburg">pburg</a>
%      <a href="matlab:help pcov">pcov</a>
%      <a href="matlab:help pmcov">pmcov</a>
%      <a href="matlab:help pyulear">pyulear</a>
%      <a href="matlab:help periodogram">periodogram</a>
%      <a href="matlab:help pmtm">pmtm</a>
%      <a href="matlab:help pwelch">pwelch</a>
%
%   Hopts = PSDOPTS(Hs) returns a PSD options object in Hopts for the PSD
%   estimator specified in Hs.
%
%   Valid PSD estimators are:
%           <a href="matlab:help spectrum.periodogram">periodogram</a>    <a href="matlab:help spectrum.mcov">mcov</a>
%           <a href="matlab:help spectrum.welch">welch</a>          <a href="matlab:help spectrum.mtm">mtm</a>
%           <a href="matlab:help spectrum.burg">burg</a>           <a href="matlab:help spectrum.yulear">yulear</a>
%           <a href="matlab:help spectrum.cov">cov</a>
%
%   Hopts contains the following properties:
%
%   Property            Valid values and description
%   ---------           ----------------------------
%   FreqPoints          [ {All} | User Defined ]
%                       Full implies full nyquist range and dynamically
%                       creates NFFT property. User Defined dynamically
%                       creates FrequencyVector property and allows the
%                       user to specify frequencies to evaluate the psd at.
%
%   NFFT                [ Auto | {Nextpow2} -or- a positive integer ]
%                       Number of FFT points. Auto uses the maximum of 256
%                       or the input (or segment for Welch) length.
%                       Nextpow2 is the same as Auto, but uses the next
%                       power of 2. This property only applies when you set
%                       the FreqPoints property to 'All'.
%
%   FrequencyVector     [ vector of real numeric doubles less than Fs ]
%                       Specify a vector of frequencies at which to
%                       evaluate the psd. This property only applies when
%                       you set the FreqPoints property to 'User Defined'.
%
%   NormalizedFrequency [ {true} | false ]
%                       False indicates that the frequency units are in
%                       Hertz.
%
%   Fs                  [ {Normalized} -or- a positive double ]
%                       Sampling frequency which can be specified only when
%                       'NormalizedFrequency' is set to false. 
%
%   SpectrumType        [ {Onesided} | Twosided ]  
%                       Onesided indicates that the total signal power is
%                       contained in half the Nyquist range.
%
%   CenterDC            [ {false} | true ]
%                       True shifts the data and frequency values so that
%                       the DC component is at the center of the spectrum.  
%                       Setting CenterDC to true forces the SpectrumType
%                       property to 'TwoSided'.
%
%   ConfLevel           [ a value between 0 and 1 ]
%                       Specifies the confidence level for computing the
%                       confidence interval, which is an estimate of the
%                       error in the calculated PSD. The confidence level is
%                       between 0 and 1. A value of 0.95 corresponds to 95%
%                       confidence level.
%
%   ConfInterval        [ N x 2 matrix ]
%                       Confidence interval. N is either NFFT or the length
%                       of the FrequencyVector. The confidence interval is
%                       usually set by the psd method of PSD estimators.
%                       The first column contains the lower bound while the
%                       second column contains the upper bound. The
%                       confidence interval contains values with a specific
%                       probability of being the actual psd. For example a
%                       95% confidence interval corresponds to the region
%                       around the nominal curve that have a 95%
%                       probability of being the true psd.
%
%   The Hopts object can be passed in as an input argument to the <a href="matlab:help spectrum/psd">psd</a> 
%   method.
%
%   Hopts = PSDOPTS(Hs,X) uses the data specified in X to return data
%   specific default options in Hopts.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help for the PSDOPTS method.

% [EOF]
