function h = lpstop(varargin)
%LPSTOP   Construct a LPSTOP object.
%   H = LPSTOP(N,Fstop,Astop,Fs) constructs a lowpass filter specifications
%   object with stopband-edge specs.
%
%   N is the filter order and must be a positive integer.
%
%   Fstop is the stopband-edge frequency and must be a positive scalar
%   between 0 and 1 if no sampling frequency is specified or between 0 and
%   Fs/2 if a sampling frequency Fs is specified.
%
%   Astop is the minimum stopband attenuation and it must be a positive
%   scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.


%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(0,4,nargin,'struct'));
    
h = fspecs.lpstop;

constructor(h,varargin{:});

h.ResponseType = 'Lowpass with stopband-edge specifications';



% [EOF]
