function h = ahpcutoff(varargin)
%AHPCUTOFF   Construct an AHPCUTOFF object.
%   H = AHPCUTOFF(N,Wc) Constructs an analog highpass filter design
%   specifications object H.
%
%   N is the filter order, and must be a positive integer.
%
%   Wc is the cutoff frequency, in radians-per-second.


%   Author(s): R. Losada
%   Copyright 1988-2003 The MathWorks, Inc.

error(nargchk(0,2,nargin,'struct'));

h = fspecs.ahpcutoff;

constructor(h,varargin{:});

h.ResponseType = 'Analog highpass with cutoff';

% [EOF]
