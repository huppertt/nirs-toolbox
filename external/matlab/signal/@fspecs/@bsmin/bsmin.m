function h = bsmin(varargin)
%BSMIN   Construct a BSMIN object.
%   H = BPMIN(Fpass1,Fstop1,Fstop2,Fpass2,Apass1,Astop,Apass2,Fs)
%   constructs a minimum-order bandstop filter specifications object.
%
%   Fpass1 is the lower passband-edge frequency and must be a positive
%   scalar   between 0 and 1 if no sampling frequency is specified or
%   between 0 and Fs/2 if a sampling frequency Fs is specified.
%
%   Fstop1 is the lower stopband-edge frequency and must be a positive
%   scalar greater than Fpass1 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Fstop2 is the higher stopband-edge frequency and must be a positive
%   scalar greater than Fstop1 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Fpass2 is the higher passband-edge frequency and must be a positive
%   scalar greater than Fstop2 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Apass1 is the maximum lower-passband deviation in dB. It must be a
%   positive scalar.
%
%   Astop is the minimum stopband attenuation in dB. It must be a positive
%   scalar.
%
%   Apass2 is the maximum higher-passband deviation in dB. It must be a
%   positive scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

h = fspecs.bsmin;
respstr = 'Minimum-order bandstop';
fstart = 1;
fstop = 4;
nargsnoFs = 7;
fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});


% [EOF]
