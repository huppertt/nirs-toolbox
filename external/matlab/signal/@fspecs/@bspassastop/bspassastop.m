function h = bspassastop(varargin)
%BSPASSASTOP   Construct a BSPASSASTOP object.
%   H = BSPASSASTOP(N,Fpass1,Fpass2,Apass,Astop,Fs) constructs a bandstop
%   filter specifications object with passband-edge specifications and
%   stopband attenuation.
%
%   N is the filter order and must be an even positive integer.
%
%   Fpass1 is the lower passband-edge frequency and must be a positive
%   scalar between 0 and 1 if no sampling frequency is specified or between
%   0 and Fs/2 if a sampling frequency Fs is specified.
%
%   Fpass2 is the higher passband-edge frequency and must be a positive
%   scalar larger than Fpass1 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
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

h = fspecs.bspassastop;
respstr = 'Bandstop with passband-edge specifications and stopband attenuation.';
fstart = 2;
fstop = 3;
nargsnoFs = 5;
fsconstructor(h,respstr,fstart,fstop,nargsnoFs,varargin{:});

% [EOF]
