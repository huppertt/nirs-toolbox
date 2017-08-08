function this = pseudospectrum(varargin)
%PSEUDOSPECTRUM   Pseudo Spectrum.
%
%   DSPDATA.PSEUDOSPECTRUM is not recommended.
%   Use the following functions instead:
%      <a href="matlab:help pmusic">pmusic</a> 
%      <a href="matlab:help rootmusic">rootmusic</a> 
%      <a href="matlab:help peig">peig</a>
%      <a href="matlab:help rooteig">rooteig</a> 
%
%   H = DSPDATA.PSEUDOSPECTRUM(DATA) instantiates an object H with its data
%   property set to DATA. DATA represents power and therefore must contain
%   real and positive values. DATA can be a vector or a matrix where each
%   column represents an independent trial. A corresponding frequency
%   vector is automatically generated in the range of [0, pi]. Fs defaults
%   to "Normalized".
%
%   H = DSPDATA.PSEUDOSPECTRUM(DATA,FREQUENCIES) sets the frequency vector
%   to FREQUENCIES in the object.  The length of the vector FREQUENCIES
%   must equal the length of the columns of DATA.
%
%   H = DSPDATA.PSEUDOSPECTRUM(...,'Fs',Fs) sets the sampling frequency to
%   Fs.  If FREQUENCIES is not specified, the frequency vector defaults to
%   [0, Fs/2].  See the NOTE below for more details.
%
%   H = DSPDATA.PSEUDOSPECTRUM(...,'SpectrumRange',SPECTRUMRANGE) sets the
%   SpectrumRange property to the string specified by SPECTRUMRANGE, which
%   can be either 'half' or 'whole'.
%
%   H = DSPDATA.PSEUDOSPECTRUM(...,'CenterDC',true) indicates that the
%   Data's DC value is centered in the vector. Setting CenterDC to true
%   automatically sets the 'SpectrumRange' to 'whole'. 
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
%   EXAMPLE: Use eigenanalysis to estimate the pseudospectrum of a noisy
%            % sinusoidal signal with two frequency components. Then store
%            % the results in a PSEUDOSPECTRUM data object and plot it.
%
%            Fs = 32e3;   t  = 0:1/Fs:2.96; 
%            x = cos(2*pi*t*1.24e3) + cos(2*pi*t*10e3) + randn(size(t)); 
%            P = pmusic(x,4);
%            hps = dspdata.pseudospectrum(P,'Fs',Fs); % Create the data object.
%            plot(hps);                               % Plot the pseudospectrum.

%   Author(s): P. Pacheco
%   Copyright 1988-2012 The MathWorks, Inc.

narginchk(0,8);

% Create and initialize object.
this = dspdata.pseudospectrum;
set(this,'Name','Pseudospectrum');

% Construct a metadata object.
set(this,'Metadata',dspdata.powermetadata);
set(this.Metadata,'FrequencyUnits','Hz');
set(this.Metadata,'DataUnits','volts^2');

% Initialize Data and Frequencies with defaults or user specified values.
initialize(this,varargin{:});

% [EOF]
