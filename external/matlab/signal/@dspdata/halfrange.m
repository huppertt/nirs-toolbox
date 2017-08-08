function halfrange(this)
%HALFRANGE   Spectrum calculated over half the Nyquist interval.
%
%   HALFRANGE is not recommended.
%   Use the 'half' option in the following functions instead:
%      <a href="matlab:help peig">peig</a>
%      <a href="matlab:help pmusic">pmusic</a>
%
%   HALFRANGE(H) converts the power spectrum specified by the data object
%   (<a href="matlab:help dspdata">dspdata</a>) H to a spectrum calculated over half the Nyquist interval.
%   The relevant properties such as, Frequencies and SpectrumRange, are
%   updated to reflect the new frequency range.
%
%   NOTE: No check is made to ensure that the data is symmetric, i.e., it
%   is assumed that the spectrum is from a real signal, and therefore only
%   half the data points are used.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help for HALFRANGE method.

% [EOF]
