function h = bpmin(varargin)
%BPMIN   Construct a BPMIN object.
%   H = BPMIN(Fstop1,Fpass1,Fpass2,Fstop2,Astop1,Apass,Astop2,Fs)
%   constructs a minimum-order bandpass filter specifications object.
%
%   Fstop1 is the lower stopband-edge frequency and must be a positive
%   scalar   between 0 and 1 if no sampling frequency is specified or
%   between 0 and Fs/2 if a sampling frequency Fs is specified.
%
%   Fpass1 is the lower passband-edge frequency and must be a positive
%   scalar greater than Fstop1 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Fpass2 is the higher passband-edge frequency and must be a positive
%   scalar greater than Fpass1 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Fstop2 is the higher stopband-edge frequency and must be a positive
%   scalar greater than Fpass2 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Astop1 is the minimum lower-stopband attenuation in dB. It must be a
%   positive scalar.
%
%   Apass is the maximum passband deviation in dB. It must be a positive
%   scalar.
%
%   Astop2 is the minimum higher-stopband attenuation in dB. It must be a
%   positive scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.bpmin;
respstr = 'Minimum-order bandpass';
fstart = 1;
fstop = 4;
nargsnoFs = 7;
fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});



% [EOF]
