function this = hilbmin(varargin)
%HILBMIN   Construct a HILBMIN object.

%   Author(s): P. Costa
%   Copyright 2005 The MathWorks, Inc.

this = fspecs.hilbmin;

this.ResponseType = 'Minimum-order Hilbert Transformer';

this.setspecs(varargin{:});


% [EOF]
