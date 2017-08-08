function this = psgaussnsym(varargin)
%PSRCOSNSYM Construct a PSRCOSNSYM object

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.psgaussnsym;

this.ResponseType = 'Gaussian pulse shaping with filter length in symbols';

this.NumberOfSymbols = 6;

this.setspecs(varargin{:});

% [EOF]
