function [Pxx,w,msg,units,Sxx,options,msgobj] = arspectra(method,x,p,varargin)
%ARSPECTRA Power Spectral Density estimate via a specified parametric method.
%   Pxx = ARSPECTRA(METHOD,X,P) returns the order P parametric PSD estimate
%   of the columns of X in the columns of Pxx.  It uses the method specified
%   in METHOD.  METHOD can be one of: 'arcov', 'arburg', 'armcov' and
%   'aryule'.  If X is a vector, it will be converted to a column vector
%   before processing.
%
%   For real signals, ARSPECTRA returns the one-sided PSD by default; for 
%   complex signals, it returns the two-sided PSD.  Note that a one-sided PSD
%   contains the total power of the input signal.
%
%   Pxx = ARSPECTRA(...,NFFT) specifies the FFT length used to calculate the
%   PSD estimates.  For real X, Pxx has length (NFFT/2+1) if NFFT is even, 
%   and (NFFT+1)/2 if NFFT is odd.  For complex X, Pxx always has length NFFT.
%   If empty, the default NFFT is 256.
%
%   [Pxx,W] = ARSPECTRA(...) returns the vector of normalized angular
%   frequencies, W, at which the PSD is estimated.  W has units of 
%   radians/sample.  For real signals, W spans the interval [0,Pi] when NFFT
%   is even and [0,Pi) when NFFT is odd.  For complex signals, W always 
%   spans the interval [0,2*Pi).
%
%   [Pxx,F] = ARSPECTRA(...,Fs) specifies a sampling frequency Fs in Hz and
%   returns the power spectral density in units of power per Hz.  F is a
%   vector of frequencies, in Hz, at which the PSD is estimated.  For real 
%   signals, F spans the interval [0,Fs/2] when NFFT is even and [0,Fs/2)
%   when NFFT is odd.  For complex signals, F always spans the interval 
%   [0,Fs).  If Fs is empty, [], the sampling frequency defaults to 1 Hz.  
%
%   [Pxx,W] = ARSPECTRA(...,'twosided') returns the PSD over the interval
%   [0,2*Pi), and [Pxx,F] = PBURG(...,Fs,'twosided') returns the PSD over
%   the interval [0,Fs).  Note that 'onesided' may be optionally specified,
%   but is only valid for real X.  The optional, trailing, string argument
%   'twosided' or 'onesided' may be placed in any position in the input 
%   argument list after the input argument P.
%
%   [...,MSG] = ARSPECTRA(...) returns an optional output string argument, MSG,
%   which contains the error message string, if any.
%
%   [...,MSG,UNITS] = ARSPECTRA(...) returns an optional output string argument, 
%   UNITS, which contains the string describing the units of the frequency 
%   vector returned.
%
%   [...,Sxx] = ARSPECTRA(...) returns the Power Spectrum (PS) estimate Sxx
%   via the specified parametric method.
%
%   ARSPECTRA complies with the general specs. for all AR PSD functions.

%   Copyright 1988-2014 The MathWorks, Inc.

nargs = 3; % Number of explicit args (not included in varargin)
narginchk(nargs,8);

if any(strcmp(varargin, 'whole'))
    warning(message('signal:arspectra:InvalidRange', '''whole''', '''twosided'''));
elseif any(strcmp(varargin, 'half'))
    warning(message('signal:arspectra:InvalidRange', '''half''', '''onesided'''));
end

% Generate an options structure with all the info. we need
if nargin == nargs,
   varargin = {}; 
end
isreal_x = isreal(x);
[options,msg,msgobj] = arspectra_options(isreal_x,varargin{:}); 

% Cast to enforce precision rules
options.nfft = double(options.nfft);
options.Fs = double(options.Fs);

if ~isempty(msg),
   Pxx = [];
   Sxx = [];
   w = [];
   units = '';
   return
end

% Determine the AR model
[a,v] = feval(method,x,p); 

% Compute the power spectrum via freqz. Always use the 'whole' option in order
% to be able to return the Nyquist point.
if isempty(options.Fs)
    Fs = {};
else
    Fs = {options.Fs};
end

if isvector(x)
    [h,w] = freqz(1,a,options.nfft(:),'whole',Fs{:});
else
    % initialize h, w
    [h,w] = freqz(1,a(1,:),options.nfft(:),'whole',Fs{:});

    % pre-initialize size of h
    h(:,2:size(a,1)) = 0;

    % fill in a spectrum for each channel
    for i=2:size(a,1)
      h(:,i) = freqz(1,a(i,:),options.nfft(:),'whole',Fs{:});
    end
end  
    
if any(v<0)
    error(message('signal:arspectra:MustBePositive'));
end

% This is the power spectrum [Power] (input variance*power of AR model)
Sxx = bsxfun(@times,v,abs(h).^2); 

% Compute the 1-sided or 2-sided PSD [Power/freq], or mean-square [Power].
% Also, compute the corresponding frequency and frequency units.
[Pxx,w,units] = computepsd(Sxx,w,options.range,options.nfft,options.Fs,'psd');

if isa(Pxx,'single')
  % Cast to enforce precision rules.
  w = single(w);
end  

%-----------------------------------------------------------------------------------
function [options,msg,msgobj] = arspectra_options(isreal_x,varargin)
%ARSPECTRA_OPTIONS   Parse the optional inputs to the ARSPECTRA function.
%   ARSPECTRA_OPTIONS returns a structure, OPTIONS, with following fields:
%
% Inputs:
%   isreal_x       - flag which is set to 1 if x is real and 0 if x is complex
%   varagin        - contains optional inputs to ARSPECTRA, such as NFFT, Fs, 
%                    and the string 'onesided' or 'twosided'
%
% Outputs:
%   options.nfft   - number of freq. points at which the psd is estimated
%   options.Fs     - sampling freq. if any
%   options.range  - 'onesided' or 'twosided' psd

% Generate defaults   
options.nfft = 256;
options.Fs = []; % Work in radians/sample
if nargin > 1 && isnumeric(varargin{1})
    nfft_freqs = varargin{1}; 
else
    nfft_freqs = options.nfft;
end
if isreal_x && (length(nfft_freqs) <= 1),
   options.range = 'onesided';
else
   options.range = 'twosided';
end

[options,msg,msgobj] = psdoptions(isreal_x,options,varargin{:});

if (length(nfft_freqs)>1 && strcmpi(options.range,'onesided'))
    warning(message('signal:arspectra:InconsistentRangeOption'));
    options.range = 'twosided';
end

% [EOF] arspectra.m
