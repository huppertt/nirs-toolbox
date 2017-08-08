function this = pssqrtrcosnsym(varargin)
%PSRCOSNSYM Construct a PSRCOSNSYM object

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.pssqrtrcosnsym;

this.ResponseType = 'Square root raised cosine pulse shaping with filter length in symbols';

this.NumberOfSymbols = 6;

this.setspecs(varargin{:});

% [EOF]
