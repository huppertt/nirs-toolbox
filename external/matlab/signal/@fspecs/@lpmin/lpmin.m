function h = lpmin(varargin)
%LPMIN   Construct a LPMIN object.
%   H = LPMIN(Fpass,Fstop,Apass,Astop,Fs) constructs a minimum-order
%   lowpass filter specifications object.
%
%   Fpass is the passband-edge frequency and must be a positive scalar
%   between 0 and 1 if no sampling frequency is specified or between 0 and
%   Fs/2 if a sampling frequency Fs is specified.
%
%   Fstop is the stopband-edge frequency and must be a positive scalar
%   greater than Fpass and  between 0 and 1 if no sampling frequency is
%   specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Apass is the maximum passband deviation in dB. It must be a positive
%   scalar.
%
%   Astop is the minimum stopband attenuation in dB. It must be a positive
%   scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.


h = fspecs.lpmin;
respstr = 'Minimum-order lowpass';
fstart = 1;
fstop = 2;
nargsnoFs = 4;
fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});



% [EOF]
