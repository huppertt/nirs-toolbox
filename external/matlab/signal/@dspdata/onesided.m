function onesided(this)
%ONESIDED   Spectrum calculated over half the Nyquist interval.
%
%   ONESIDED is not recommended.
%   Use the 'onesided' option in the following functions instead:
%      <a href="matlab:help periodogram">periodogram</a>
%      <a href="matlab:help pwelch">pwelch</a>
%      <a href="matlab:help pburg">pburg</a>
%      <a href="matlab:help pcov">pcov</a>
%      <a href="matlab:help pmcov">pmcov</a>
%      <a href="matlab:help pyulear">pyulear</a>
%      <a href="matlab:help pmtm">pmtm</a>
%
%   ONESIDED(H) converts the power spectrum specified by the data object
%   (<a href="matlab:help dspdata">dspdata</a>) H to a spectrum calculated over half the Nyquist interval
%   containing the full power. The relevant properties such as, Frequencies
%   and SpectrumType, are updated to reflect the new frequency range.
%
%   NOTE: No check is made to ensure that the data is symmetric, i.e., it
%   is assumed that the spectrum is from a real signal, and therefore only
%   half the data points are used.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help for the ONESIDED method.

% [EOF]
