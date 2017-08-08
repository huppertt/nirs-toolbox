function plot(this)
%PLOT   Plot DSP data (DSPDATA) objects.
%
%   DSPDATA.PLOT is not recommended.  
%   Use the following functions instead:
%      <a href="matlab:help periodogram">periodogram</a>
%      <a href="matlab:help pwelch">pwelch</a>
%      <a href="matlab:help pburg">pburg</a>
%      <a href="matlab:help pcov">pcov</a>
%      <a href="matlab:help pmcov">pmcov</a>
%      <a href="matlab:help pyulear">pyulear</a>
%      <a href="matlab:help pmtm">pmtm</a>
%      <a href="matlab:help peig">peig</a>
%      <a href="matlab:help pmusic">pmusic</a>
%
%   PLOT(H) plots the data in the <a href="matlab:help dspdata">dspdata</a> object H.
%
%   EXAMPLE: Use the periodogram to estimate the power spectral density of
%            % a noisy sinusoidal signal with two frequency components.
%            % Then store the results in PSD data object and plot it.
%
%            Fs = 32e3;   t = 0:1/Fs:2.96;
%            x = cos(2*pi*t*1.24e3)+ cos(2*pi*t*10e3)+ randn(size(t));
%            Pxx = periodogram(x);
%            hpsd = dspdata.psd(Pxx,'Fs',Fs); % Create a PSD data object.
%            plot(hpsd);                      % Plot the PSD.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

% Help for dspdata PLOT method.

% [EOF]
