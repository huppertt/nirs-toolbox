function this = msspectrum(varargin)
%MSSPECTRUM   Mean-square spectrum.
%
%   DSPDATA.MSSPECTRUM is not recommended.
%   Use <a href="matlab:help periodogram">periodogram</a> and <a href="matlab:help pwelch">pwelch</a> instead.
%
%   H = DSPDATA.MSSPECTRUM(DATA) instantiates an object H with its data
%   property set to the mean-square spectrum specified in DATA. DATA
%   represents power and therefore must contain real and positive values.
%   DATA can be a vector or a matrix where each column represents an
%   independent trial. A corresponding frequency vector is automatically
%   generated in the range of [0, pi]. Fs defaults to "Normalized".
%
%   The mean-squared spectrum is intended for discrete spectra. Unlike the
%   power spectral density (PSD), the peaks in the mean-square spectrum
%   reflect the power in the signal at a given frequency.
%
%   H = DSPDATA.MSSPECTRUM(DATA,FREQUENCIES) sets the frequency vector to
%   FREQUENCIES in the data object returned in H.  The length of the vector
%   FREQUENCIES must equal the length of the columns of DATA.
%
%   H = DSPDATA.MSSPECTRUM(...,'Fs',Fs) sets the sampling frequency to Fs.
%   If FREQUENCIES is not specified the frequency vector defaults to
%   [0,Fs/2]. See the NOTE below for more details.
%
%   H = DSPDATA.MSSPECTRUM(...,'SpectrumType',SPECTRUMTYPE) sets the
%   SpectrumType property to the string specified by SPECTRUMTYPE, which
%   can be either 'onesided' or 'twosided'.
%
%   H = DSPDATA.MSSPECTRUM(...,'CenterDC',true) indicates that the Data's
%   DC value is centered in the vector. Setting CenterDC to true
%   automatically sets the 'SpectrumType' to 'twosided'. 
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
%   EXAMPLE: Use FFT to calculate mean-square spectrum of a noisy
%            % sinusoidal signal with two frequency components. Then store
%            % the results in an MSSPECTRUM data object and plot it.
%
%            Fs = 32e3;   t = 0:1/Fs:2.96;
%            x = cos(2*pi*t*1.24e3) + cos(2*pi*t*10e3) + randn(size(t));
%            X = fft(x);
%            P = (abs(X)/length(x)).^2;    % Compute the mean-square.
%        
%            hms = dspdata.msspectrum(P,'Fs',Fs,'SpectrumType','twosided'); 
%            plot(hms);                    % Plot the mean-square spectrum.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

narginchk(0,12);  

% Create object and set the properties specific to this object.
this = dspdata.msspectrum;
set(this,'Name','Mean-Square Spectrum');

% Construct a metadata object.
set(this,'Metadata',dspdata.powermetadata);
set(this.Metadata,'FrequencyUnits','Hz');
set(this.Metadata,'DataUnits','volts^2');

% Initialize Data and Frequencies with defaults or user specified values.
initialize(this,varargin{:});

% [EOF]
