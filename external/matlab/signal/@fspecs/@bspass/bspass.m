function h = bspass(varargin)
%BSPASS   Construct a BSPASS object.
%   H = BSPASS(N,Fpass1,Fpass2,Apass,Fs) constructs a bandstop filter
%   specifications object with passband-edge specifications.
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
%   Apass is the maximum passband deviation and it must be a positive
%   scalar.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(0,5,nargin,'struct'));
h = fspecs.bspass;
constructor(h,varargin{:});
h.ResponseType = 'Bandstop with passband-edge specifications.';

% [EOF]
