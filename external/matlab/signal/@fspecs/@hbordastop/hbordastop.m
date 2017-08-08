function this = hbordastop(varargin)
%HBORDASTOP   Construct a HBORDASTOP object.

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

this = fspecs.hbordastop;

this.ResponseType = 'Halfband with filter order and stopband attenuation';

this.setspecs(varargin{:});

% [EOF]
