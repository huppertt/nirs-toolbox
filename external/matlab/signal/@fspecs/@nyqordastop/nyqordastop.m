function this = nyqordastop(varargin)
%NYQORDASTOP   Construct a NYQORDASTOP object.
%   NYQORDASTOP 

%   Author(s): J. Schickler
%   Copyright 1988-2003 The MathWorks, Inc.

this = fspecs.nyqordastop;

this.ResponseType = 'Nyquist with filter order and stopband attenuation';

this.setspecs(varargin{:});

% [EOF]
