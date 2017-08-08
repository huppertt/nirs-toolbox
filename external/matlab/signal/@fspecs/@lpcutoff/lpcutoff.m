function h = lpcutoff(varargin)
%LPCUTOFF   Construct a LPCUTOFF object.
%   H = LPCUTOFF(N,Fcutoff,Fs) constructs a lowpass filter specifications
%   object with a cutoff frequency.
%
%   N is the filter order and must be a positive integer.
%
%   Fcutoff is the cutoff frequency and must be a positive scalar between 0
%   and 1 if no sampling frequency is specified or between 0 and Fs/2 if a
%   sampling frequency Fs is specified.
%
%   Fs is the sampling frequency. If Fs is not specified, normalized
%   frequency is assumed. If Fs is specified, it must be a positive scalar.

%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(0,3,nargin,'struct'));


h = fspecs.lpcutoff;

constructor(h,varargin{:});

h.ResponseType = 'Lowpass with cutoff';

% [EOF]
