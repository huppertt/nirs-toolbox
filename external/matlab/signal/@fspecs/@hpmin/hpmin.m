function h = hpmin(varargin)
%HPMIN   Construct a HPMIN object.
%   H = HPMIN(Fstop,Fpass,Astop,Apass,Fs) constructs a minimum-order
%   highpass filter specifications object.
%
%   Fstop is the stopband-edge frequency and must be a positive scalar
%   between 0 and 1 if no sampling frequency is specified or between 0 and
%   Fs/2 if a sampling frequency Fs is specified.
%
%   Fpass is the passband-edge frequency and must be a positive scalar
%   greater than Fstop and between 0 and 1 if no sampling frequency is
%   specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Astop is the minimum stopband attenuation in dB. It must be a positive
%   scalar.
%
%   Apass is the maximum passband deviation in dB. It must be a positive
%   scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

% Override factory defaults inherited from lowpass
if nargin < 1,
    varargin{1} = .45;
end
if nargin < 2,
    varargin{2} = .55;
end


h = fspecs.hpmin;
respstr = 'Minimum-order highpass';
fstart = 1;
fstop = 2;
nargsnoFs = 4;
fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});


% [EOF]
