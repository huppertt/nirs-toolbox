function hopts = msspectrumopts(this,x)
%MSSPECTRUMOPTS   Mean-square spectrum (MSS) options object.
%
%   MSSPECTRUMOPTS is not recommended.  Use <a href="matlab:help periodogram">periodogram</a> or <a href="matlab:help pwelch">pwelch</a> instead.
%
%   Hopts = MSSPECTRUMOPTS(Hs) returns an MSS options object in Hopts for
%   the MSS estimator specified in Hs.
%
%    Valid MSS estimators:
%           <a href="matlab:help spectrum.periodogram">periodogram</a>
%           <a href="matlab:help spectrum.welch">welch</a>
%
%   Hopts contains the following properties:
%
%   Property            Valid values and description
%   ---------           ----------------------------
%   FreqPoints          [ {All} | User Defined ]
%                       Full implies full nyquist range and dynamically
%                       creates NFFT property. User Defined dynamically
%                       creates FrequencyVector property and allows the
%                       user to specify frequencies to evaluate the
%                       mean-square spectrum at.
%
%   NFFT                [ Auto | {Nextpow2} -or- a positive integer ]
%                       Number of FFT points. Auto uses the maximum of 256
%                       or the input (or segment for Welch) length.
%                       Nextpow2 is the same as Auto, but uses the next
%                       power of 2.
%
%   FrequencyVector     [ vector of real numeric doubles less than Fs ]
%                       Specify a vector of frequencies at which to
%                       evaluate the mean-square spectrum.
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
%   The Hopts object can be passed in as an input argument to the <a href="matlab:help spectrum/msspectrum">msspectrum</a> 
%   method.
%
%   Hopts = MSSPECTRUMOPTS(Hs,X) uses the data specified in X to return
%   data specific default options in Hopts.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help for the MSSPECTRUMOPTS method.

% [EOF]
