function h = bpstop(varargin)
%BPSTOP   Construct a BPSTOP object.
%   H = BPSTOP(N,Fstop1,Fstop2,Astop,Fs) constructs a bandpass filter
%   specifications object with stopband-edge specifications.
%
%   N is the filter order and must be an even positive integer.
%
%   Fstop1 is the lower stopband-edge frequency and must be a positive
%   scalar between 0 and 1 if no sampling frequency is specified or between
%   0 and Fs/2 if a sampling frequency Fs is specified.
%
%   Fstop2 is the higher stopband-edge frequency and must be a positive
%   scalar larger than Fstop1 and between 0 and 1 if no sampling frequency
%   is specified or between 0 and Fs/2 if a sampling frequency Fs is
%   specified.
%
%   Astop is the minimum stopband attenuation and it must be a positive
%   scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(0,5,nargin,'struct'));

h = fspecs.bpstop;
constructor(h,varargin{:});
h.ResponseType = 'Bandpass with stopband-edge specifications.';


% [EOF]
