function this = psd(varargin)
%PSD   Power Spectral Density (PSD).
%
%   DSPDATA.PSD is not recommended.  Use the following functions instead:
%      <a href="matlab:help periodogram">periodogram</a>
%      <a href="matlab:help pwelch">pwelch</a>
%      <a href="matlab:help pburg">pburg</a>
%      <a href="matlab:help pcov">pcov</a>
%      <a href="matlab:help pmcov">pmcov</a>
%      <a href="matlab:help pyulear">pyulear</a>
%      <a href="matlab:help pmtm">pmtm</a>
%
%   H = DSPDATA.PSD(DATA) instantiates a data object H with its data
%   property set to DATA. DATA represents power and therefore must contain
%   real and positive values. DATA can be a vector or a matrix where each
%   column represents an independent trial. A corresponding frequency
%   vector is automatically generated in the range of [0, pi]. Fs
%   defaults to "Normalized".
%
%   The power spectral density is intended for continuous spectra. Note
%   that unlike the mean-squared spectrum (MSS), in this case the peaks in
%   the spectra do not reflect the power at a given frequency. Instead,
%   the integral of the PSD over a given frequency band computes the
%   average power in the signal over such frequency band. See the help on
%   AVGPOWER for more information.
%
%   H = DSPDATA.PSD(DATA,FREQUENCIES) sets the frequency vector to
%   FREQUENCIES in the data object returned in H.  The length of the vector
%   FREQUENCIES must equal the length of the columns of DATA.
%
%   H = DSPDATA.PSD(...,'Fs',Fs) sets the sampling frequency to Fs.  If
%   FREQUENCIES is not specified the frequency vector defaults to [0,Fs/2].
%   See the NOTE below for more details.
%
%   H = DSPDATA.PSD(...,'SpectrumType',SPECTRUMTYPE) sets the SpectrumType
%   property to the string specified by SPECTRUMTYPE, which can be either
%   'onesided' or 'twosided'.
%
%   H = DSPDATA.PSD(...,'CenterDC',true) indicates that the Data's DC value
%   is centered in the vector. Setting CenterDC to true automatically sets
%   the 'SpectrumType' to 'twosided'.
%
%   If no frequency vector is specified the default frequency vector is
%   generated according to the setting of 'CenterDC'.  If a frequency
%   vector is specified then 'CenterDC' should be set to match the
%   frequency vector (and data) specified.  To modify this property use the
%   <a href="matlab:help dspdata/centerdc">centerdc</a> method.
%
%   NOTE: If the spectrum data specified was calculated over "half" the
%   Nyquist interval and you don't specify a corresponding frequency
%   vector, then the default frequency vector will assume that the number
%   of points in the "whole" FFT was even.  Also, the plot option to
%   convert to a "whole" spectrum will assume the original "whole" FFT
%   length was even.
%
%   EXAMPLE: Use the periodogram to estimate the power spectral density of
%            % a noisy sinusoidal signal with two frequency components. Then
%            % store the results in PSD data object and plot it.
%
%            Fs = 32e3;   t = 0:1/Fs:2.96;
%            x = cos(2*pi*t*1.24e3)+ cos(2*pi*t*10e3)+ randn(size(t));
%            Pxx = periodogram(x);
%            hpsd = dspdata.psd(Pxx,'Fs',Fs); % Create a PSD data object.
%            plot(hpsd);                      % Plot the PSD.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

narginchk(0,12);

% Create object and set the properties specific to this object.
this = dspdata.psd;
set(this,'Name','Power Spectral Density');

% Construct a metadata object.
set(this,'Metadata',dspdata.powermetadata);
set(this.Metadata,...
    'FrequencyUnits','Hz',...
    'DataUnits','volts^2/Hz');

% Initialize Data and Frequencies with defaults or user specified values.
initialize(this,varargin{:});

% [EOF]
