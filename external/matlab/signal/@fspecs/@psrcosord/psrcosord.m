function this = psrcosord(varargin)
%PSRCOSORD Construct a PSRCOSORD object

%   Copyright 2008 The MathWorks, Inc.

this = fspecs.psrcosord;

this.ResponseType = 'Raised cosine pulse shaping with filter order';

this.FilterOrder = 48;  % (SamplesPerSymbol * NumberOfSymbols = 8*6)

this.setspecs(varargin{:});

% [EOF]
