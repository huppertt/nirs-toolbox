function this = hbordntw(varargin)
%HBORDNTW   Construct a HBORDNTW object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

this = fspecs.hbordntw;

this.ResponseType = 'Halfband with filter order and transition width';

this.setspecs(varargin{:});

% [EOF]
