function this = hbord(varargin)
%HBORD   Construct a HBORD object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

this = fspecs.hbord;

this.ResponseType = 'Halfband with filter order';

this.setspecs(varargin{:});

% [EOF]
