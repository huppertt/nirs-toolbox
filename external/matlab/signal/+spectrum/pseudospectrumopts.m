function hopts = pseudospectrumopts(this,x)
%PSEUDOSPECTRUMOPTS   Pseudospectrum options object.
%
%   PSEUDOSPECTRUMOPTS is not recommended.  
%   Use the following functions instead:
%      <a href="matlab:help pmusic">pmusic</a> 
%      <a href="matlab:help rootmusic">rootmusic</a> 
%      <a href="matlab:help peig">peig</a>
%      <a href="matlab:help rooteig">rooteig</a> 
%
%   Hopts = PSEUDOSPECTRUMOPTS(Hs) returns a pseudospectrum options object
%   in Hopts for the pseudospectrum estimator specified in Hs.
%
%    Valid pseudospectrum estimators are:
%           <a href="matlab:help spectrum.eigenvector">eigenvector</a>
%           <a href="matlab:help spectrum.music">music</a>
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
%                       pseudospectrum at.
%
%   NFFT                [ Auto | {Nextpow2} -or- a positive integer ]
%                       Number of FFT points. Auto uses the maximum of 256
%                       or the input (or segment for Welch) length.
%                       Nextpow2 is the same as Auto, but uses the next
%                       power of 2.
%
%   FrequencyVector     [ vector of real numeric doubles less than Fs ]
%                       Specify a vector of frequencies at which to
%                       evaluate the pseudospectrum.
%
%   NormalizedFrequency [ {true} | false ]
%                       False indicates that the frequency units are in
%                       Hertz.
%
%   Fs                  [ {Normalized} -or- a positive double ]
%                       Sampling frequency which can be specified only when
%                       'NormalizedFrequency' is set to false. 
%
%   SpectrumRange        [ {Half} | Whole ]  
%                       Half indicates that only half the Nyquist range is
%                       used.
%
%   The Hopts object can be passed in as an input argument to the 
%   <a href="matlab:help spectrum/pseudospectrum">pseudospectrum</a> method.
%
%   Hopts = PSEUDOSPECTRUMOPTS(Hs,X) uses the data specified in X to return
%   data specific default options in Hopts.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help for the PSEUDOSPECTRUMOPTS method.

% [EOF]
